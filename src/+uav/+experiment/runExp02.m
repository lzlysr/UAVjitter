function result = runExp02(cfg)
%RUNEXP02 Run the paired six-method jitter sweep.

experimentTimer = tic;
[cfg, calibration] = uav.experiment.prepare(cfg);
experimentIndex = 2;
jitterPoints = cfg.experiments.jitterSweepDeg;
pointCount = numel(jitterPoints);
trialCount = cfg.monteCarlo.nTrials;
periodCount = cfg.experiments.periodCounts.jitterSweep;
methodCount = numel(cfg.experiments.methodNames);
raw = allocateRaw(trialCount, periodCount, pointCount, methodCount);
seeds = cell(1, pointCount);
for pointIndex = 1:pointCount
    pointCfg = cfg;
    pointCfg.attitude.stdDeg = jitterPoints(pointIndex) * ones(size(cfg.attitude.stdDeg));
    pointCfg = uav.config.validate(pointCfg);
    pointCfg.beamforming.epsilonS = cfg.beamforming.epsilonS;
    [pointRaw, seeds{pointIndex}] = uav.experiment.collectMethods( ...
        pointCfg, experimentIndex, pointIndex, periodCount);
    raw = assignPoint(raw, pointRaw, pointIndex);
end
sinrSummary = uav.experiment.summarizePoints( ...
    raw.sinrLinear, cfg, experimentIndex);
gainSummary = uav.experiment.summarizeFinitePoints(raw.desiredGain, cfg);
leakageSummary = uav.experiment.summarizeFinitePoints( ...
    raw.interferenceLeakage, cfg);
summaryTable = sweepTable(jitterPoints, cfg.experiments.methodNames, ...
    sinrSummary, gainSummary, leakageSummary);

fig = uav.plot.newFigure(cfg);
layout = tiledlayout(fig, 'flow');
nexttile(layout); hold on;
for methodIndex = 1:methodCount
    lower = squeeze(sinrSummary.ci95Db(1, :, methodIndex)).';
    upper = squeeze(sinrSummary.ci95Db(2, :, methodIndex)).';
    center = sinrSummary.medianDb(:, methodIndex);
    errorbar(jitterPoints, center, center - lower, upper - center, '-o', ...
        'LineWidth', cfg.plot.lineWidth, 'MarkerSize', cfg.plot.markerSize);
end
xlabel('Attitude jitter RMS (deg)'); ylabel('Median output SINR (dB)'); grid on;
legend(cellstr(cfg.experiments.methodNames), 'Location', 'best');
nexttile(layout);
plot(jitterPoints, sinrSummary.lowDb, '-o', ...
    'LineWidth', cfg.plot.lineWidth, 'MarkerSize', cfg.plot.markerSize);
xlabel('Attitude jitter RMS (deg)'); ylabel('5th-percentile output SINR (dB)'); grid on;
legend(cellstr(cfg.experiments.methodNames), 'Location', 'best');
title(layout, 'Paired jitter sweep');

result = struct('experimentName', "exp02_jitter_sweep", ...
    'calibration', calibration, 'methodNames', cfg.experiments.methodNames, ...
    'seedScheme', {seeds}, 'raw', raw, ...
    'summary', struct('sinr', sinrSummary, 'desiredGain', gainSummary, ...
        'interferenceLeakage', leakageSummary));
actualCounts = struct('trials', trialCount, 'periodsPerPoint', periodCount, ...
    'points', pointCount, 'methods', methodCount);
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
raw.snapshotSignatures = zeros(trials, periods, points, methods, 4);
end

function raw = assignPoint(raw, pointRaw, pointIndex)
fields = {'sinrLinear', 'desiredGain', 'interferenceLeakage', 'runtimeSec', ...
    'incmNmse', 'incmMinimumEigenvalue', 'incmConditionNumber'};
for index = 1:numel(fields)
    raw.(fields{index})(:, :, pointIndex, :) = pointRaw.(fields{index});
end
raw.snapshotSignatures(:, :, pointIndex, :, :) = pointRaw.snapshotSignatures;
end

function output = sweepTable(points, names, sinr, gain, leakage)
[pointGrid, methodGrid] = ndgrid(points, 1:numel(names));
output = table(pointGrid(:), reshape(names(methodGrid), [], 1), ...
    sinr.medianDb(:), reshape(sinr.ci95Db(1, :, :), [], 1), ...
    reshape(sinr.ci95Db(2, :, :), [], 1), sinr.lowDb(:), ...
    gain.median(:), leakage.median(:), ...
    'VariableNames', {'jitterRmsDeg', 'method', 'medianSinrDb', ...
    'ci95LowDb', 'ci95HighDb', 'p05SinrDb', 'medianDesiredGain', ...
    'medianInterferenceLeakage'});
end
