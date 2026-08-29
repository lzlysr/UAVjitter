function result = runClassicJitterSweep(cfg)
%RUNCLASSICJITTERSWEEP 经典方法在不同姿态抖动RMS下的输出SINR退化实验。
%   该实验独立于固定五实验主线，方法组合完全由
%   cfg.experiments.classicJitterMethodNames决定，不自动加入run_all。

experimentTimer = tic;
[cfg, calibration] = uav.experiment.prepare(cfg);
experimentIndex = 6;
jitterPoints = cfg.experiments.jitterSweepDeg;
methodNames = reshape(cfg.experiments.classicJitterMethodNames, 1, []);
methodLabels = classicMethodLabels(methodNames);
pointCount = numel(jitterPoints);
trialCount = cfg.monteCarlo.nTrials;
% 独立实验沿用原抖动扫描的周期数，避免引入仅供该实验使用的新字段。
periodCount = cfg.experiments.periodCounts.jitterSweep;
methodCount = numel(methodNames);
raw = allocateRaw(trialCount, periodCount, pointCount, methodCount);
seeds = cell(1, pointCount);

for pointIndex = 1:pointCount
    pointCfg = cfg;
    pointCfg.attitude.stdDeg = jitterPoints(pointIndex) * ...
        ones(size(cfg.attitude.stdDeg));
    pointCfg = uav.config.validate(pointCfg);
    % WCO半径由实验开始前的固定设计分布标定，扫描时不使用当前真值重算。
    pointCfg.beamforming.epsilonS = cfg.beamforming.epsilonS;
    [pointRaw, seeds{pointIndex}] = uav.experiment.collectClassicMethods( ...
        pointCfg, experimentIndex, pointIndex, periodCount, methodNames);
    raw = assignPoint(raw, pointRaw, pointIndex);
end

sinrSummary = uav.experiment.summarizePoints( ...
    raw.sinrLinear, cfg, experimentIndex);
gainSummary = uav.experiment.summarizeFinitePoints(raw.desiredGain, cfg);
leakageSummary = uav.experiment.summarizeFinitePoints( ...
    raw.interferenceLeakage, cfg);
summaryTable = classicSweepTable(jitterPoints, methodNames, ...
    sinrSummary, gainSummary, leakageSummary);

fig = uav.plot.newFigure(cfg);
layout = tiledlayout(fig, 'flow');
nexttile(layout); hold on;
% 含误差棒的中位数 95%置信区间：反映中位数估计的不确定程度，不是原始SINR样本的分布范围。
% for methodIndex = 1:methodCount
%     lower = squeeze(sinrSummary.ci95Db(1, :, methodIndex)).';
%     upper = squeeze(sinrSummary.ci95Db(2, :, methodIndex)).';
%     center = sinrSummary.medianDb(:, methodIndex);
%     errorbar(jitterPoints, center, center - lower, upper - center, '-o', ...
%         'LineWidth', cfg.plot.lineWidth, 'MarkerSize', cfg.plot.markerSize);
% end
% 只含中位数
for methodIndex = 1:methodCount
    center = sinrSummary.medianDb(:, methodIndex);
    [lineStyle, marker] = classicMethodStyle(methodNames(methodIndex));
    curve = plot(jitterPoints, center, ...
        'LineStyle', lineStyle, ...
        'Marker', marker, ...
        'LineWidth', cfg.plot.lineWidth, ...
        'MarkerSize', cfg.plot.markerSize-4);
    if methodNames(methodIndex) == "oracle"
        % Oracle是理想性能上界，使用黑色虚线与其他比较方法区分。
        set(curve, 'Color', [0, 0, 0], 'LineStyle', '--');
    end
end
xlabel('姿态抖动RMS（°）', 'FontSize', cfg.plot.fontSize);
ylabel('输出SINR（dB）', 'FontSize', cfg.plot.fontSize); grid on;
set(gca, 'FontSize', cfg.plot.fontSize, 'LineWidth', cfg.plot.lineWidth-1);
box on;
ylim([-30, 30]);
legend(cellstr(methodLabels), 'Location', 'bestoutside', ...
    'FontSize', cfg.plot.fontSize);
% nexttile(layout);
% plot(jitterPoints, sinrSummary.lowDb, '-o', ...
%     'LineWidth', cfg.plot.lineWidth, 'MarkerSize', cfg.plot.markerSize);
% xlabel('姿态抖动RMS（°）', 'FontSize', cfg.plot.fontSize);
% ylabel('输出SINR第5百分位（dB）', 'FontSize', cfg.plot.fontSize); grid on;
% set(gca, 'FontSize', cfg.plot.fontSize, 'LineWidth', cfg.plot.lineWidth-1);
% box on;
% legend(cellstr(methodLabels), 'Location', 'bestoutside', ...
%     'FontSize', cfg.plot.fontSize);
% title(layout, '经典方法在姿态抖动下的输出SINR退化');

result = struct( ...
    'experimentName', "exp06_classic_method_jitter_sweep", ...
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
end

function [lineStyle, marker] = classicMethodStyle(name)
%CLASSICMETHODSTYLE 为每种经典方法返回固定且互不重复的绘图样式。
catalogNames = ["oracle", "nsb", "frozen_mvdr", "scm_mvdr", "nlms", ...
    "rls", "pi", "dl_mvdr", "esp", "wco", "static_incm_mvdr"];
lineStyles = ["--", "-", "--", "-.", "-.", "-", "--", "-.", ...
    ":", "-", "--"];
markers = ["none", "o", "s", "^", "v", "d", "x", "+", "p", ...
    "h", "*"];
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
    "rls", "pi", "dl_mvdr", "esp", "wco", "static_incm_mvdr"];
catalogLabels = ["理想上界", "NSB", "冻结初始MVDR", "MVDR", ...
    "NLMS", "RLS", "PI", "DL", "ESP", "WCO", ...
    "INCM"];
[found, locations] = ismember(names, catalogNames);
if ~all(found)
    error('uav:experiment:InvalidClassicMethodSelection', ...
        'A classic method has no display label.');
end
labels = catalogLabels(locations);
end

function output = classicSweepTable(points, names, sinr, gain, leakage)
[pointGrid, methodGrid] = ndgrid(points, 1:numel(names));
output = table(pointGrid(:), reshape(names(methodGrid), [], 1), ...
    sinr.medianDb(:), reshape(sinr.ci95Db(1, :, :), [], 1), ...
    reshape(sinr.ci95Db(2, :, :), [], 1), sinr.lowDb(:), ...
    gain.median(:), leakage.median(:), ...
    'VariableNames', {'jitterRmsDeg', 'method', 'medianSinrDb', ...
    'ci95LowDb', 'ci95HighDb', 'p05SinrDb', 'medianDesiredGain', ...
    'medianInterferenceLeakage'});
end
