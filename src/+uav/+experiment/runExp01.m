function result = runExp01(cfg)
%RUNEXP01 实验1：验证姿态—方向—阵列流形动态模型。
%   本实验不比较六种波束形成算法，而是回答模型是否按预期工作：
%   1) 三轴姿态偏差是否随时间变化；
%   2) 世界固定方向在阵列坐标系中是否产生时变方向余弦；
%   3) 名义流形与真实流形的失配和阵列响应如何变化；
%   4) 独立零抖动场景是否严格保持不变。

experimentTimer = tic;
[cfg, calibration] = uav.experiment.prepare(cfg);
experimentIndex = 1;
periodCount = cfg.experiments.periodCounts.modelValidation;
scenarioSubstream = uav.experiment.substream(cfg, experimentIndex, 1, 1, 0, "scenario");
scenario = uav.scenario.generate(cfg, ...
    uav.util.makeStream(cfg.random.baseSeed, scenarioSubstream));
sourceCount = cfg.derived.J + 1;
mismatch = zeros(sourceCount, periodCount);
arrayResponse = zeros(sourceCount, periodCount);
for periodIndex = 1:periodCount
    steering = scenario.steeringTruth(:, :, periodIndex);
    nominal = cfg.derived.nominalSteering;
    % 除以sqrt(M)，得到不随阵元总数线性放大的归一化导向失配。
    mismatch(:, periodIndex) = vecnorm(steering - nominal, 2, 1).' / sqrt(cfg.derived.M);
    % |a_nom^H a_true|^2/M^2：名义波束观察真实方向时的归一化响应。
    arrayResponse(:, periodIndex) = ...
        abs(sum(conj(nominal) .* steering, 1)).'.^2 / cfg.derived.M^2;
end

zeroCfg = cfg;
% 零抖动验证使用独立场景子流，避免直接复用已有数组而形成“自证”。
zeroCfg.attitude.stdDeg = zeros(size(cfg.attitude.stdDeg));
zeroCfg = uav.config.validate(zeroCfg);
zeroCfg.beamforming.epsilonS = cfg.beamforming.epsilonS;
zeroSubstream = uav.experiment.substream(cfg, experimentIndex, 2, 1, 0, "scenario");
zeroScenario = uav.scenario.generate(zeroCfg, ...
    uav.util.makeStream(cfg.random.baseSeed, zeroSubstream));
zeroMismatch = max(abs(zeroScenario.steeringTruth(:, :, 1:periodCount) - ...
    zeroCfg.derived.nominalSteering), [], 'all');
zeroDirectionVariation = max(abs(zeroScenario.directionCosines(:, :, 1:periodCount) - ...
    zeroScenario.directionCosines(:, :, 1)), [], 'all');
if zeroMismatch > cfg.experiments.oracleTolerance || ...
        zeroDirectionVariation > cfg.experiments.oracleTolerance
    error('uav:experiment:ZeroJitterInvarianceFailure', ...
        'The separate zero-jitter validation scenario was not invariant.');
end

timeMs = (0:periodCount - 1).' * cfg.timeline.dtSec * 1e3;
attitudeDeg = rad2deg(scenario.attitudeOffsetRad(:, 1:periodCount)).';
summaryTable = table(timeMs, attitudeDeg(:, 1), attitudeDeg(:, 2), ...
    attitudeDeg(:, 3), mismatch(1, :).', arrayResponse(1, :).', ...
    'VariableNames', {'timeMs', 'rollDeg', 'pitchDeg', 'yawDeg', ...
    'desiredMismatchNormalized', 'desiredArrayResponseNormalized'});
for sourceIndex = 1:sourceCount
    % 每个源分别导出xi_x和xi_y，CSV列名携带源序号。
    summaryTable.("xiX_source" + sourceIndex) = ...
        reshape(scenario.directionCosines(1, sourceIndex, 1:periodCount), [], 1);
    summaryTable.("xiY_source" + sourceIndex) = ...
        reshape(scenario.directionCosines(2, sourceIndex, 1:periodCount), [], 1);
end

fig = uav.plot.newFigure(cfg);
% 四个面板依次展示姿态、方向余弦、导向失配和名义阵列响应。
layout = tiledlayout(fig, 'flow');
nexttile(layout);
plot(timeMs, attitudeDeg, 'LineWidth', cfg.plot.lineWidth);
xlabel('Time (ms)'); ylabel('Attitude offset (deg)');
legend({'roll', 'pitch', 'yaw'}, 'Location', 'best'); grid on;
nexttile(layout);
plot(timeMs, reshape(permute(scenario.directionCosines(:, :, 1:periodCount), ...
    [3, 2, 1]), periodCount, []), 'LineWidth', cfg.plot.lineWidth);
xlabel('Time (ms)'); ylabel('Direction cosine'); grid on;
nexttile(layout);
plot(timeMs, mismatch.', 'LineWidth', cfg.plot.lineWidth);
xlabel('Time (ms)'); ylabel('Normalized steering mismatch'); grid on;
nexttile(layout);
plot(timeMs, arrayResponse.', 'LineWidth', cfg.plot.lineWidth);
xlabel('Time (ms)'); ylabel('Normalized nominal array response'); grid on;
title(layout, 'Model validation diagnostics');

result = struct( ...
    'experimentName', "exp01_model_validation", ...
    'calibration', calibration, ...
    'methodNames', strings(1, 0), ...
    'seedScheme', struct('scenarioSubstream', scenarioSubstream, ...
        'zeroJitterScenarioSubstream', zeroSubstream), ...
    'raw', struct('timeMs', timeMs, 'attitudeOffsetDeg', attitudeDeg, ...
        'directionCosines', scenario.directionCosines(:, :, 1:periodCount), ...
        'steeringMismatchNormalized', mismatch, ...
        'arrayResponseNormalized', arrayResponse, ...
        'zeroJitterDirectionCosines', ...
            zeroScenario.directionCosines(:, :, 1:periodCount)), ...
    'summary', struct('zeroJitterMaximumSteeringDifference', zeroMismatch, ...
        'zeroJitterMaximumDirectionVariation', zeroDirectionVariation, ...
        'arrayResponseDefinition', "|a_nominal^H a_true|^2/M^2"));
actualCounts = struct('trials', 1, 'periods', periodCount, ...
    'sources', sourceCount, 'zeroJitterPeriods', periodCount);
result = uav.experiment.finalize(result, summaryTable, fig, cfg, ...
    actualCounts, experimentTimer);
end
