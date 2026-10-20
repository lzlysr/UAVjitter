function result = exp07_ula_classic_snapshot_sweep(cfg)
%EXP07_ULA_CLASSIC_SNAPSHOT_SWEEP 在y轴ULA中执行用户指定的六方法快拍扫描。
%   一维角度theta=asin(u_y)，避免沿用UPA二维谱造成不可辨识谱脊。
%   static_incm_mvdr在本实验特指main中CovFit思路的一维移植；不改变
%   本分支其他经典实验的历史标识。原始Rec在L<M时明确记录不适用。

timer = tic;
settings = cfg.experiments.ulaClassicSnapshot;
cfg.array.Nx = 1;
cfg.array.Ny = settings.sensorCount;
cfg.experiments.snapshotSweep = settings.snapshotSweep;
% 对照main中CovFit锁定参数，而不是误用研究分支的无正则网格变体。
source = cfg.literature.guLeshem2012.upaDiagnosis;
fields = {'loadingRatio', 'fitRegularization', 'peakThresholdRatio', ...
    'peakNmsRadiusBins', 'maximumCandidateCount'};
for fi = 1:numel(fields)
    cfg.experiments.recUlaGeometry.(fields{fi}) = source.(fields{fi});
end
cfg.experiments.recUlaGeometry.sensorCount = settings.sensorCount;
cfg.experiments.recUlaGeometry.spacingOverLambda = cfg.array.spacingOverLambda(2);
cfg.experiments.recUlaGeometry.candidateGridSize = settings.angleGridSize;
cfg = uav.config.validate(cfg);
% 9阵元ULA的原姿态分位数标定会使球形约束不可行；按用户确认使用
% 原论文固定半径sqrt(0.3*M)，不改变真实姿态分布，也不按扫描结果调参。
cfg.beamforming.epsilonS = sqrt(settings.wcoSquaredRadiusScale * cfg.derived.M);
assert(cfg.beamforming.epsilonS < norm(cfg.derived.nominalSteeringReference), ...
    'WCO半径必须小于名义导向矢量范数。');
calibration = struct('epsilonS', cfg.beamforming.epsilonS, ...
    'scope', "用户确认的文献固定半径", ...
    'squaredRadiusScale', settings.wcoSquaredRadiusScale);
names = settings.methodNames;
baselineNames = ["oracle", "dl_mvdr", "esp", "wco"];
assert(isequal(names, [baselineNames, "gu_leshem2012_rec", "static_incm_mvdr"]), ...
    '本专项实验方法列表必须保持指定的六方法顺序。');
points = settings.snapshotSweep;
T = cfg.monteCarlo.nTrials;
K = cfg.experiments.periodCounts.snapshotSweep;
M = cfg.derived.M;
thetaRef = asin(cfg.derived.nominalDirectionCosines(2, 1));
step = pi / settings.angleGridSize;
theta = -pi/2 + step/2 + step*(0:settings.angleGridSize-1);
keep = abs(theta-thetaRef) > deg2rad(cfg.spectrum.desiredExclusionHalfWidthDeg);
scan = uav.array.steeringULA(theta(keep), M, cfg.array.spacingOverLambda(2));
areas = step*ones(1, sum(keep));
raw.sinrLinear = nan(T, K, numel(points), numel(names));
raw.runtimeSec = nan(size(raw.sinrLinear));
raw.status = zeros(size(raw.sinrLinear), 'int8'); % 1成功，0不适用。
raw.snapshotSignatures = zeros(T, K, numel(points), numel(names), 4);
raw.notApplicableReason = strings(numel(points), numel(names));
seeds = cell(1, numel(points));
for pointIndex = 1:numel(points)
    point = cfg;
    point.timeline.L = points(pointIndex);
    point = uav.config.validate(point);
    seeds{pointIndex} = struct('scenarioSubstreams', zeros(T,1), ...
        'snapshotSubstreams', zeros(T,K));
    for tr = 1:T
        sceneSeed = uav.experiment.substream(cfg, settings.experimentIndex, pointIndex, tr, 0, "scenario");
        scene = uav.scenario.generate(point, uav.util.makeStream(cfg.random.baseSeed, sceneSeed));
        seeds{pointIndex}.scenarioSubstreams(tr) = sceneSeed;
        for k = 1:K
            dataSeed = uav.experiment.substream(cfg, settings.experimentIndex, pointIndex, tr, k, "snapshot");
            seeds{pointIndex}.snapshotSubstreams(tr,k) = dataSeed;
            obs = uav.scenario.generatePeriod(point, scene, k, ...
                uav.util.makeStream(cfg.random.baseSeed, dataSeed));
            % 两种重构与四个基线共享同一次观测；仅Oracle和评价器读取真值。
            baseline = uav.experiment.runClassicMethods(obs, [], point, baselineNames);
            raw.sinrLinear(tr,k,pointIndex,1:4) = baseline.sinrLinear;
            raw.runtimeSec(tr,k,pointIndex,1:4) = baseline.runtimeSec;
            raw.status(tr,k,pointIndex,1:4) = 1;
            safe = uav.experiment.algorithmObservation(obs);
            for mi = 5:6
                if mi == 5 && points(pointIndex) < M
                    raw.notApplicableReason(pointIndex,mi) = "L<M，SCM秩亏，未经加载的原始Rec不适用";
                    continue;
                end
                methodTimer = tic;
                if mi == 5
                    Rh = uav.incm.guLeshem2012Rec(safe.Rx, scan, areas);
                else
                    Rh = uav.incm.ulaCovfit(safe, thetaRef, point, ...
                        cfg.spectrum.desiredExclusionHalfWidthDeg, false);
                end
                w = uav.beamforming.solve("mvdr", Rh, ...
                    cfg.derived.nominalSteeringReference, point);
                raw.runtimeSec(tr,k,pointIndex,mi) = toc(methodTimer);
                metrics = uav.metrics.evaluate(w, obs.truth, point);
                raw.sinrLinear(tr,k,pointIndex,mi) = metrics.sinrLinear;
                raw.status(tr,k,pointIndex,mi) = 1;
            end
            signature = uav.util.snapshotSignature(obs.X);
            raw.snapshotSignatures(tr,k,pointIndex,:,:) = repmat(reshape(signature,1,1,1,1,4),1,1,1,6,1);
            values = reshape(raw.sinrLinear(tr,k,pointIndex,:),1,[]);
            assert(all(isfinite(values(raw.status(tr,k,pointIndex,:)==1))), '有效方法不得有非有限SINR。');
            assert(all(values(1)+cfg.experiments.oracleTolerance >= values(isfinite(values))), ...
                '非Oracle方法不应超过理想上界。');
        end
    end
    fprintf('ULA实验7：%d/%d，L=%d完成\n', pointIndex, numel(points), points(pointIndex));
end
summary = uav.experiment.summarizePoints(raw.sinrLinear, cfg, settings.experimentIndex);
labels = ["理想上界", "DL", "ESP", "WCO", "原始Rec", "CovFit-INCM"];
fig = uav.plot.newFigure(cfg);
axesHandle = axes(fig);
hold(axesHandle, 'on');
styles = ["--", "-.", ":", "-", "-.", "--"];
markers = ["none", "+", "p", "h", "s", "*"];
for mi = 1:numel(names)
    curve = plot(axesHandle, points, summary.medianDb(:,mi), ...
        'LineStyle', styles(mi), 'Marker', markers(mi), ...
        'LineWidth', cfg.plot.lineWidth, 'MarkerSize', cfg.plot.markerSize);
    if mi == 1, curve.Color = [0,0,0]; end
end
xlabel(axesHandle, '每周期快拍数 L', 'FontSize', cfg.plot.fontSize);
ylabel(axesHandle, '输出SINR（dB）', 'FontSize', cfg.plot.fontSize);
set(axesHandle, 'FontSize', cfg.plot.fontSize);
grid(axesHandle, 'on'); box(axesHandle, 'on');
legend(axesHandle, cellstr(labels), 'Location', 'bestoutside', 'FontSize', cfg.plot.fontSize);
[pgrid,mgrid] = ndgrid(points,1:numel(names));
tableOut = table(pgrid(:), reshape(names(mgrid),[],1), summary.medianDb(:), ...
    summary.validCounts(:), reshape(raw.notApplicableReason,[],1), ...
    'VariableNames', {'snapshotCount','method','medianSinrDb','validCount','notApplicableReason'});
result = struct('experimentName', "exp07_ula_classic_snapshot_sweep", ...
    'methodNames', names, 'methodLabels', labels, 'scanPoints', points, ...
    'calibration', calibration, 'raw', raw, 'summary', struct('sinr',summary), ...
    'seedScheme', {seeds}, 'ulaReferenceRad', thetaRef, ...
    'ulaGridAnglesRad', theta, 'ulaIntegrationMask', keep, ...
    'covfitDefinition', "main锁定CovFit参数的一维ULA移植，非本分支旧static积分");
counts = struct('trials',T,'periodsPerPoint',K,'points',numel(points),'methods',numel(names));
result = uav.experiment.finalize(result,tableOut,fig,cfg,counts,timer);
end
