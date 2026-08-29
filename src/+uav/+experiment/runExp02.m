function result = runExp02(cfg)
%RUNEXP02 实验2：配对扫描三轴姿态抖动RMS。
%   每个扫描点只改变cfg.attitude.stdDeg，其余几何、功率、算法参数和离线
%   epsilonS保持不变。同一点内全部所选方法共享场景与快拍；方法组合由
%   cfg.experiments.exp02MethodNames给出，不影响其他实验的固定六方法合同。

experimentTimer = tic;
[cfg, calibration] = uav.experiment.prepare(cfg);
experimentIndex = 2;
jitterPoints = cfg.experiments.jitterSweepDeg;
methodNames = reshape(cfg.experiments.exp02MethodNames, 1, []);
pointCount = numel(jitterPoints);
trialCount = cfg.monteCarlo.nTrials;
periodCount = cfg.experiments.periodCounts.jitterSweep;
methodCount = numel(methodNames);
raw = allocateRaw(trialCount, periodCount, pointCount, methodCount);
seeds = cell(1, pointCount);
for pointIndex = 1:pointCount
    % 本实验按批准方案令三轴使用相同RMS扫描值。
    pointCfg = cfg;
    pointCfg.attitude.stdDeg = jitterPoints(pointIndex) * ones(size(cfg.attitude.stdDeg));
    pointCfg = uav.config.validate(pointCfg);
    % validate会重建派生量，因此显式恢复实验开始时一次标定的epsilonS。
    pointCfg.beamforming.epsilonS = cfg.beamforming.epsilonS;
    [pointRaw, seeds{pointIndex}] = uav.experiment.collectMethods( ...
        pointCfg, experimentIndex, pointIndex, periodCount, methodNames);
    raw = assignPoint(raw, pointRaw, pointIndex);
end
sinrSummary = uav.experiment.summarizePoints( ...
    raw.sinrLinear, cfg, experimentIndex);
gainSummary = uav.experiment.summarizeFinitePoints(raw.desiredGain, cfg);
leakageSummary = uav.experiment.summarizeFinitePoints( ...
    raw.interferenceLeakage, cfg);
% 失配是场景指标而不是方法指标。先在每次试验内部跨周期取中位数，
% 再统计不同独立试次的分布，避免把时间相关周期误当成独立样本。
trialMedianMismatch = median( ...
    raw.desiredSteeringMismatchNormalized, 2);
mismatchSummary = uav.experiment.summarizeFinitePoints( ...
    trialMedianMismatch, cfg);
mismatchSummary.aggregationUnit = "trial-median-over-periods";
summaryTable = sweepTable(jitterPoints, methodNames, ...
    sinrSummary, gainSummary, leakageSummary, mismatchSummary);

fig = uav.plot.newFigure(cfg);
% 三个面板依次展示中位SINR、SINR低分位和式（3-10）期望导向失配。
layout = tiledlayout(fig, 'flow');
nexttile(layout); hold on;
for methodIndex = 1:methodCount
    lower = squeeze(sinrSummary.ci95Db(1, :, methodIndex)).';
    upper = squeeze(sinrSummary.ci95Db(2, :, methodIndex)).';
    center = sinrSummary.medianDb(:, methodIndex);
    errorbar(jitterPoints, center, center - lower, upper - center, '-o', ...
        'LineWidth', cfg.plot.lineWidth, 'MarkerSize', cfg.plot.markerSize);
end
xlabel('姿态抖动RMS（°）', 'FontSize', cfg.plot.fontSize);
ylabel('输出SINR中位数（dB）', 'FontSize', cfg.plot.fontSize); grid on;
set(gca, 'FontSize', cfg.plot.fontSize, 'LineWidth', cfg.plot.lineWidth-1);
legend(cellstr(methodNames), 'Location', 'best', ...
    'FontSize', cfg.plot.fontSize);
nexttile(layout);
plot(jitterPoints, sinrSummary.lowDb, '-o', ...
    'LineWidth', cfg.plot.lineWidth, 'MarkerSize', cfg.plot.markerSize);
xlabel('姿态抖动RMS（°）', 'FontSize', cfg.plot.fontSize);
ylabel('输出SINR第5百分位（dB）', 'FontSize', cfg.plot.fontSize); grid on;
set(gca, 'FontSize', cfg.plot.fontSize, 'LineWidth', cfg.plot.lineWidth-1);
legend(cellstr(methodNames), 'Location', 'best', ...
    'FontSize', cfg.plot.fontSize);
nexttile(layout);
mismatchMedian = mismatchSummary.median(:, 1);
mismatchLow = mismatchSummary.low(:, 1);
mismatchHigh = mismatchSummary.p95(:, 1);
errorbar(jitterPoints, mismatchMedian, ...
    mismatchMedian - mismatchLow, mismatchHigh - mismatchMedian, '-o', ...
    'LineWidth', cfg.plot.lineWidth, 'MarkerSize', cfg.plot.markerSize);
xlabel('姿态抖动RMS（°）', 'FontSize', cfg.plot.fontSize);
ylabel('归一化期望导向失配程度', ...
    'FontSize', cfg.plot.fontSize); grid on;
ylim([0,1]);
set(gca, 'FontSize', cfg.plot.fontSize, 'LineWidth', cfg.plot.lineWidth-1);
% title(layout, '姿态抖动RMS配对扫描');

result = struct('experimentName', "exp02_jitter_sweep", ...
    'calibration', calibration, 'methodNames', methodNames, ...
    'seedScheme', {seeds}, 'raw', raw, ...
    'summary', struct('sinr', sinrSummary, 'desiredGain', gainSummary, ...
        'interferenceLeakage', leakageSummary, ...
        'desiredSteeringMismatch', mismatchSummary));
actualCounts = struct('trials', trialCount, 'periodsPerPoint', periodCount, ...
    'points', pointCount, 'methods', methodCount);
result = uav.experiment.finalize(result, summaryTable, fig, cfg, ...
    actualCounts, experimentTimer);
end

function raw = allocateRaw(trials, periods, points, methods)
%ALLOCATERAW 按[试次×周期×扫描点×方法]预分配，避免循环中增长数组。
fields = {'sinrLinear', 'desiredGain', 'interferenceLeakage', 'runtimeSec'};
for index = 1:numel(fields)
    raw.(fields{index}) = zeros(trials, periods, points, methods);
end
fields = {'incmNmse', 'incmMinimumEigenvalue', 'incmConditionNumber'};
for index = 1:numel(fields)
    raw.(fields{index}) = nan(trials, periods, points, methods);
end
raw.snapshotSignatures = zeros(trials, periods, points, methods, 4);
% 期望导向失配只由姿态场景决定，不沿方法维重复存储。
raw.desiredSteeringMismatchNormalized = zeros(trials, periods, points);
end

function raw = assignPoint(raw, pointRaw, pointIndex)
%ASSIGNPOINT 把单一扫描点的三维结果写入总结果的指定point切片。
fields = {'sinrLinear', 'desiredGain', 'interferenceLeakage', 'runtimeSec', ...
    'incmNmse', 'incmMinimumEigenvalue', 'incmConditionNumber'};
for index = 1:numel(fields)
    raw.(fields{index})(:, :, pointIndex, :) = pointRaw.(fields{index});
end
raw.snapshotSignatures(:, :, pointIndex, :, :) = pointRaw.snapshotSignatures;
raw.desiredSteeringMismatchNormalized(:, :, pointIndex) = ...
    pointRaw.desiredSteeringMismatchNormalized;
end

function output = sweepTable(points, names, sinr, gain, leakage, mismatch)
%SWEEPTABLE 将“扫描点×方法”的矩阵统计展开为长表，便于CSV分析。
[pointGrid, methodGrid] = ndgrid(points, 1:numel(names));
mismatchMedian = repmat(mismatch.median(:, 1), 1, numel(names));
mismatchLow = repmat(mismatch.low(:, 1), 1, numel(names));
mismatchHigh = repmat(mismatch.p95(:, 1), 1, numel(names));
output = table(pointGrid(:), reshape(names(methodGrid), [], 1), ...
    sinr.medianDb(:), reshape(sinr.ci95Db(1, :, :), [], 1), ...
    reshape(sinr.ci95Db(2, :, :), [], 1), sinr.lowDb(:), ...
    gain.median(:), leakage.median(:), mismatchMedian(:), ...
    mismatchLow(:), mismatchHigh(:), ...
    'VariableNames', {'jitterRmsDeg', 'method', 'medianSinrDb', ...
    'ci95LowDb', 'ci95HighDb', 'p05SinrDb', 'medianDesiredGain', ...
    'medianInterferenceLeakage', ...
    'medianDesiredSteeringMismatchNormalized', ...
    'p05DesiredSteeringMismatchNormalized', ...
    'p95DesiredSteeringMismatchNormalized'});
end
