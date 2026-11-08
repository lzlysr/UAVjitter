function result = runClassicSweep(cfg, scanKind)
%RUNCLASSICSWEEP 实验6/7/8共用的配对扫描、统计和导出流程。
%   jitter只改变三轴相同RMS；snapshots只改变L并固定配置的三轴RMS。
%   snr只改变期望输入SNR；L和三轴RMS由实验8专用配置固定，不改INR。
%   每个扫描点内全部方法共享观测，不自动加入run_all。

experimentTimer = tic;
cfg = uav.config.validate(cfg);
switch string(scanKind)
    case "jitter"
        experimentIndex = 6;
        scanPoints = cfg.experiments.jitterSweepDeg;
        methodNames = reshape(cfg.experiments.classicJitterMethodNames, 1, []);
        periodCount = cfg.experiments.periodCounts.jitterSweep;
        xLabel = '姿态抖动RMS（°）';
        pointVariable = 'jitterRmsDeg';
        experimentName = "exp06_classic_method_jitter_sweep";
    case "snapshots"
        experimentIndex = 7;
        scanPoints = cfg.experiments.snapshotSweep;
        methodNames = reshape(cfg.experiments.classicSnapshotMethodNames, 1, []);
        periodCount = cfg.experiments.periodCounts.snapshotSweep;
        xLabel = '每周期快拍数 L';
        pointVariable = 'snapshotCount';
        experimentName = "exp07_classic_method_snapshot_sweep";
    case "snr"
        experimentIndex = 8;
        scanPoints = cfg.experiments.snrSweepDb;
        methodNames = reshape(cfg.experiments.classicSnrMethodNames, 1, []);
        periodCount = cfg.experiments.periodCounts.snrSweep;
        xLabel = '输入SNR（dB）';
        pointVariable = 'inputSnrDb';
        experimentName = "exp08_classic_method_snr_sweep";
        % 在标定及生成观测之前设置固定条件，最终MAT也保存实际使用值。
        cfg.timeline.L = cfg.experiments.classicSnrSnapshotCount;
        cfg.attitude.stdDeg = cfg.experiments.classicSnrJitterStdDeg;
        cfg = uav.config.validate(cfg);
    otherwise
        error('uav:experiment:InvalidClassicScan', '扫描类型必须为jitter、snapshots或snr。');
end
if cfg.derived.arrayType == "ULA" && ~any(methodNames=="wco")
    % 经典INCM都接标准MVDR；未选择WCO时不应被其不可行半径阻断。
    % 不改变配置半径；若选择WCO，则由统一标定入口使用显式ULA半径配置。
    calibration=struct('epsilonS',cfg.beamforming.epsilonS, ...
        'scope',"未选择WCO，标准MVDR对照不需要鲁棒半径标定");
else
    [cfg, calibration] = uav.experiment.prepare(cfg);
end
methodLabels = classicMethodLabels(methodNames);
if cfg.derived.arrayType == "ULA"
    overlap=cfg.derived.ulaGeometry.interferersInsideDesiredSector;
    if any(overlap)
        warning('uav:experiment:UlaSectorOverlap', ...
            'ULA名义干扰%s落入期望排除区；保留main参数，不自动改变分区或源方向。',mat2str(find(overlap)));
    end
end
pointCount = numel(scanPoints);
trialCount = cfg.monteCarlo.nTrials;
methodCount = numel(methodNames);
raw = allocateRaw(trialCount, periodCount, pointCount, methodCount);
seeds = cell(1, pointCount);

for pointIndex = 1:pointCount
    pointCfg = cfg;
    switch string(scanKind)
        case "jitter"
            pointCfg.attitude.stdDeg = scanPoints(pointIndex) * ...
                ones(size(cfg.attitude.stdDeg));
        case "snapshots"
            pointCfg.timeline.L = scanPoints(pointIndex);
        case "snr"
            pointCfg.signal.snrDb = scanPoints(pointIndex);
    end
    pointCfg = uav.config.validate(pointCfg);
    % WCO半径由实验开始前的固定设计分布标定，扫描时不使用当前真值重算。
    pointCfg.beamforming.epsilonS = cfg.beamforming.epsilonS;
    [pointRaw, seeds{pointIndex}] = uav.experiment.collectClassicMethods( ...
        pointCfg, experimentIndex, pointIndex, periodCount, methodNames);
    raw = assignPoint(raw, pointRaw, pointIndex);
    fprintf('%s：扫描点 %d/%d 完成（%g）\n', experimentName, ...
        pointIndex, pointCount, scanPoints(pointIndex));
end

sinrSummary = uav.experiment.summarizePoints( ...
    raw.sinrLinear, cfg, experimentIndex);
gainSummary = uav.experiment.summarizeFinitePoints(raw.desiredGain, cfg);
leakageSummary = uav.experiment.summarizeFinitePoints( ...
    raw.interferenceLeakage, cfg);
summaryTable = classicSweepTable(scanPoints, methodNames, ...
    sinrSummary, gainSummary, leakageSummary, pointVariable);
summaryTable.validCount = sinrSummary.validCounts(:);
summaryTable.notApplicableCount = reshape(sum(sum(raw.methodStatus=="not-applicable",1),2),[],1);
reasons = strings(pointCount,methodCount);
for p=1:pointCount
    for m=1:methodCount
        entries=unique(reshape(raw.notApplicableReason(:,:,p,m),[],1));
        reasons(p,m)=strjoin(entries(strlength(entries)>0),"；");
    end
end
summaryTable.notApplicableReason = reasons(:);

fig = uav.plot.newFigure(cfg);
layout = tiledlayout(fig, 'flow');
nexttile(layout); hold on;
% 图中只画中位数；95%中位数置信区间仍保存到MAT和CSV中。
for methodIndex = 1:methodCount
    center = sinrSummary.medianDb(:, methodIndex);
    [lineStyle, marker] = classicMethodStyle(methodNames(methodIndex));
    curve = plot(scanPoints, center, ...
        'LineStyle', lineStyle, ...
        'Marker', marker, ...
        'LineWidth', cfg.plot.lineWidth, ...
        'MarkerSize', cfg.plot.markerSize-4);
    if methodNames(methodIndex) == "oracle"
        % Oracle是理想性能上界，使用黑色虚线与其他比较方法区分。
        set(curve, 'Color', [0, 0, 0], 'LineStyle', '--');
    end
end
xlabel(xLabel, 'FontSize', cfg.plot.fontSize);
ylabel('输出SINR（dB）', 'FontSize', cfg.plot.fontSize); grid on;
set(gca, 'FontSize', cfg.plot.fontSize, 'LineWidth', cfg.plot.lineWidth-1);
box on;
% 不固定纵轴范围，避免较差方法被裁掉而掩盖实际性能。
legend(cellstr(methodLabels), 'Location', 'bestoutside', ...
    'FontSize', cfg.plot.fontSize, 'Interpreter','none');

result = struct( ...
    'experimentName', experimentName, ...
    'scanKind', string(scanKind), 'scanPoints', scanPoints, ...
    'calibration', calibration, ...
    'methodNames', methodNames, ...
    'methodLabels', methodLabels, ...
    'seedScheme', {seeds}, ...
    'raw', raw, ...
    'summary', struct('sinr', sinrSummary, ...
        'desiredGain', gainSummary, ...
        'interferenceLeakage', leakageSummary));
actualCounts = struct('trials', trialCount, ...
    'periodsPerPoint', periodCount, 'points', pointCount, ...
    'methods', methodCount);
result = uav.experiment.finalize(result, summaryTable, fig, cfg, ...
    actualCounts, experimentTimer);
end

function raw = allocateRaw(trials, periods, points, methods)
fields = {'sinrLinear', 'desiredGain', 'interferenceLeakage', 'runtimeSec'};
for index = 1:numel(fields)
    raw.(fields{index}) = zeros(trials, periods, points, methods);
end
fields = {'incmNmse', 'incmMinimumEigenvalue', 'incmConditionNumber'};
for index = 1:numel(fields)
    raw.(fields{index}) = nan(trials, periods, points, methods);
end
raw.desiredSteeringMismatchNormalized = zeros(trials, periods, points);
raw.snapshotSignatures = zeros(trials, periods, points, methods, 4);
raw.methodStatus = strings(trials,periods,points,methods);
raw.notApplicableReason = strings(trials,periods,points,methods);
end

function raw = assignPoint(raw, pointRaw, pointIndex)
fields = {'sinrLinear', 'desiredGain', 'interferenceLeakage', 'runtimeSec', ...
    'incmNmse', 'incmMinimumEigenvalue', 'incmConditionNumber'};
for index = 1:numel(fields)
    raw.(fields{index})(:, :, pointIndex, :) = pointRaw.(fields{index});
end
raw.desiredSteeringMismatchNormalized(:, :, pointIndex) = ...
    pointRaw.desiredSteeringMismatchNormalized;
raw.snapshotSignatures(:, :, pointIndex, :, :) = pointRaw.snapshotSignatures;
raw.methodStatus(:,:,pointIndex,:) = pointRaw.methodStatus;
raw.notApplicableReason(:,:,pointIndex,:) = pointRaw.notApplicableReason;
end

function [lineStyle, marker] = classicMethodStyle(name)
%CLASSICMETHODSTYLE 为每种经典方法返回固定且互不重复的绘图样式。
catalogNames = ["oracle", "nsb", "frozen_mvdr", "scm_mvdr", "nlms", ...
    "rls", "pi", "dl_mvdr", "esp", "wco", "incm_covfit", "incm_integral_processed", "incm_covfit_refined", "incm_integral_origin", ...
    "incm_capon_point","incm_urglq"];
lineStyles = ["--", "-", "--", "-.", "-.", "-", "--", "-.", ...
    ":", "-", "--", "-.", ":", "--", ":", "-."];
markers = ["none", "o", "s", "^", "v", "d", "x", "+", "p", ...
    "h", "*", "s", ">", "<", "o", "d"];
[found, location] = ismember(string(name), catalogNames);
if ~found
    error('uav:experiment:InvalidClassicMethodSelection', ...
        'A classic method has no plotting style.');
end
lineStyle = lineStyles(location);
marker = markers(location);
end

function labels = classicMethodLabels(names)
catalogNames = ["oracle", "nsb", "frozen_mvdr", "scm_mvdr", "nlms", ...
    "rls", "pi", "dl_mvdr", "esp", "wco", "incm_covfit", "incm_integral_processed", "incm_covfit_refined", "incm_integral_origin", ...
    "incm_capon_point","incm_urglq"];
catalogLabels = ["理想上界", "NSB", "冻结初始MVDR", "MVDR", ...
    "NLMS", "RLS", "PI", "DL", "ESP", "WCO", ...
    "incm_covfit", "incm_integral_processed", "incm_covfit_refined", "incm_integral_origin", ...
    "incm_capon_point","incm_urglq"];
[found, locations] = ismember(names, catalogNames);
if ~all(found)
    error('uav:experiment:InvalidClassicMethodSelection', ...
        'A classic method has no display label.');
end
labels = catalogLabels(locations);
end

function output = classicSweepTable(points, names, sinr, gain, leakage, pointVariable)
[pointGrid, methodGrid] = ndgrid(points, 1:numel(names));
output = table(pointGrid(:), reshape(names(methodGrid), [], 1), ...
    sinr.medianDb(:), reshape(sinr.ci95Db(1, :, :), [], 1), ...
    reshape(sinr.ci95Db(2, :, :), [], 1), sinr.lowDb(:), ...
    gain.median(:), leakage.median(:), ...
    'VariableNames', {pointVariable, 'method', 'medianSinrDb', ...
    'ci95LowDb', 'ci95HighDb', 'p05SinrDb', 'medianDesiredGain', ...
    'medianInterferenceLeakage'});
end
