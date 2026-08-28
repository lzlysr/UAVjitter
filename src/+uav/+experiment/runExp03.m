function result = runExp03(cfg)
%RUNEXP03 Run the paired snapshot-count sweep.

experimentTimer = tic;
[cfg, calibration] = uav.experiment.prepare(cfg);
experimentIndex = 3;
snapshotPoints = cfg.experiments.snapshotSweep;
pointCount = numel(snapshotPoints);
trialCount = cfg.monteCarlo.nTrials;
periodCount = cfg.experiments.periodCounts.snapshotSweep;
methodCount = numel(cfg.experiments.methodNames);
raw = allocateRaw(trialCount, periodCount, pointCount, methodCount);
seeds = cell(1, pointCount);
for pointIndex = 1:pointCount
    pointCfg = cfg;
    pointCfg.timeline.L = snapshotPoints(pointIndex);
    pointCfg = uav.config.validate(pointCfg);
    pointCfg.beamforming.epsilonS = cfg.beamforming.epsilonS;
    [pointRaw, seeds{pointIndex}] = uav.experiment.collectMethods( ...
        pointCfg, experimentIndex, pointIndex, periodCount);
    raw = assignPoint(raw, pointRaw, pointIndex);
end
sinrSummary = uav.experiment.summarizePoints(raw.sinrLinear, cfg, experimentIndex);
nmseSummary = uav.experiment.summarizeFinitePoints(raw.incmNmse, cfg);
minimumEigenvalueSummary = uav.experiment.summarizeFinitePoints( ...
    raw.incmMinimumEigenvalue, cfg);
conditionNumberSummary = uav.experiment.summarizeFinitePoints( ...
    raw.incmConditionNumber, cfg);
runtimeSummary = uav.experiment.summarizeFinitePoints(raw.runtimeSec, cfg);
summaryTable = snapshotTable(snapshotPoints, cfg.experiments.methodNames, ...
    sinrSummary, nmseSummary, minimumEigenvalueSummary, ...
    conditionNumberSummary, runtimeSummary);

fig = uav.plot.newFigure(cfg);
layout = tiledlayout(fig, 'flow');
nexttile(layout);
semilogx(snapshotPoints, sinrSummary.medianDb, '-o', ...
    'LineWidth', cfg.plot.lineWidth, 'MarkerSize', cfg.plot.markerSize);
xlabel('Snapshots per period'); ylabel('Median output SINR (dB)'); grid on;
legend(cellstr(cfg.experiments.methodNames), 'Location', 'best');
nexttile(layout);
semilogx(snapshotPoints, nmseSummary.median, '-o', ...
    'LineWidth', cfg.plot.lineWidth, 'MarkerSize', cfg.plot.markerSize);
xlabel('Snapshots per period'); ylabel('Median INCM NMSE'); grid on;
legend(cellstr(cfg.experiments.methodNames), 'Location', 'best');
nexttile(layout);
semilogx(snapshotPoints, runtimeSummary.median, '-o', ...
    'LineWidth', cfg.plot.lineWidth, 'MarkerSize', cfg.plot.markerSize);
xlabel('Snapshots per period'); ylabel('Median runtime (s)'); grid on;
nexttile(layout);
semilogx(snapshotPoints, minimumEigenvalueSummary.median, '-o', ...
    'LineWidth', cfg.plot.lineWidth, 'MarkerSize', cfg.plot.markerSize);
xlabel('Snapshots per period'); ylabel('Median minimum eigenvalue'); grid on;
legend(cellstr(cfg.experiments.methodNames), 'Location', 'best');
nexttile(layout);
semilogx(snapshotPoints, conditionNumberSummary.median, '-o', ...
    'LineWidth', cfg.plot.lineWidth, 'MarkerSize', cfg.plot.markerSize);
xlabel('Snapshots per period'); ylabel('Median INCM condition number'); grid on;
legend(cellstr(cfg.experiments.methodNames), 'Location', 'best');
title(layout, 'Snapshot-count sweep');

result = struct('experimentName', "exp03_snapshot_sweep", ...
    'calibration', calibration, 'methodNames', cfg.experiments.methodNames, ...
    'seedScheme', {seeds}, 'raw', raw, ...
    'summary', struct('sinr', sinrSummary, 'incmNmse', nmseSummary, ...
        'incmMinimumEigenvalue', minimumEigenvalueSummary, ...
        'incmConditionNumber', conditionNumberSummary, 'runtime', runtimeSummary));
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

function output = snapshotTable(points, names, sinr, nmse, minimumEigenvalue, ...
        conditionNumber, runtime)
[pointGrid, methodGrid] = ndgrid(points, 1:numel(names));
output = table(pointGrid(:), reshape(names(methodGrid), [], 1), ...
    sinr.medianDb(:), sinr.lowDb(:), nmse.median(:), ...
    minimumEigenvalue.median(:), conditionNumber.median(:), runtime.median(:), ...
    'VariableNames', {'snapshotCount', 'method', 'medianSinrDb', 'p05SinrDb', ...
    'medianIncmNmse', 'medianIncmMinimumEigenvalue', ...
    'medianIncmConditionNumber', 'medianRuntimeSec'});
end
