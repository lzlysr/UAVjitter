function result = runGuLeshem2012Reproduction(cfg)
%RUNGULESHEM2012REPRODUCTION 运行Gu–Leshem 2012例1的两组复现实验。
%   所有方法在每个点、每次试验共享同一组实际DOA和快拍；先在线性域
%   汇总蒙特卡洛均值，再转换为dB，以避免直接平均dB产生统计偏差。

timerValue = tic;
% 入口重新校验，防止用户在config()返回后修改字段却忘记显式调用校验器。
cfg = uav.config.validate(cfg);
validateConfiguration(cfg);
paper = cfg.literature.guLeshem2012;
methods = reshape(string(paper.methodNames), 1, []);
trialCount = paper.monteCarloTrials;
snrPoints = reshape(paper.snrSweepDb, 1, []);
snapshotPoints = reshape(paper.snapshotSweep, 1, []);
methodCount = numel(methods);
experimentIndex = paper.experimentIndex;

raw = struct( ...
    'snrSinrLinear', zeros(trialCount, numel(snrPoints), methodCount), ...
    'snapshotSinrLinear', zeros(trialCount, numel(snapshotPoints), methodCount), ...
    'snrIncmNmse', nan(trialCount, numel(snrPoints), methodCount), ...
    'snapshotIncmNmse', nan(trialCount, numel(snapshotPoints), methodCount), ...
    'snrRuntimeSec', zeros(trialCount, numel(snrPoints), methodCount), ...
    'snapshotRuntimeSec', zeros(trialCount, numel(snapshotPoints), methodCount));

for pointIndex = 1:numel(snrPoints)
    for trialIndex = 1:trialCount
        substream = uav.experiment.substream(cfg, experimentIndex, ...
            pointIndex, trialIndex, 0, "scenario");
        obs = uav.scenario.generateGuLeshem2012Trial(cfg, ...
            snrPoints(pointIndex), paper.fixedSnapshotCount, ...
            uav.util.makeStream(cfg.random.baseSeed, substream));
        record = uav.experiment.runGuLeshem2012Methods(obs, cfg, methods);
        raw.snrSinrLinear(trialIndex, pointIndex, :) = ...
            reshape(record.sinrLinear, 1, 1, []);
        raw.snrIncmNmse(trialIndex, pointIndex, :) = ...
            reshape(record.incmNmse, 1, 1, []);
        raw.snrRuntimeSec(trialIndex, pointIndex, :) = ...
            reshape(record.runtimeSec, 1, 1, []);
    end
end

snapshotPointOffset = numel(snrPoints);
for pointIndex = 1:numel(snapshotPoints)
    for trialIndex = 1:trialCount
        substream = uav.experiment.substream(cfg, experimentIndex, ...
            snapshotPointOffset + pointIndex, trialIndex, 0, "scenario");
        obs = uav.scenario.generateGuLeshem2012Trial(cfg, ...
            paper.fixedSnrDb, snapshotPoints(pointIndex), ...
            uav.util.makeStream(cfg.random.baseSeed, substream));
        record = uav.experiment.runGuLeshem2012Methods(obs, cfg, methods);
        raw.snapshotSinrLinear(trialIndex, pointIndex, :) = ...
            reshape(record.sinrLinear, 1, 1, []);
        raw.snapshotIncmNmse(trialIndex, pointIndex, :) = ...
            reshape(record.incmNmse, 1, 1, []);
        raw.snapshotRuntimeSec(trialIndex, pointIndex, :) = ...
            reshape(record.runtimeSec, 1, 1, []);
    end
end

snrSummary = summarizeSweep(raw.snrSinrLinear, cfg, ...
    experimentIndex, 0);
snapshotSummary = summarizeSweep(raw.snapshotSinrLinear, cfg, ...
    experimentIndex, snapshotPointOffset);
summary = struct('snr', snrSummary, 'snapshot', snapshotSummary);
fig = plotReproduction(snrPoints, snapshotPoints, methods, summary, cfg);
summaryTable = reproductionTable(snrPoints, snapshotPoints, methods, summary);

result = struct( ...
    'experimentName', "exp07_gu_leshem2012_reproduction", ...
    'methodNames', methods, ...
    'paperParameters', paper, ...
    'raw', raw, ...
    'summary', summary);
actualCounts = struct( ...
    'trials', trialCount, ...
    'snrPoints', numel(snrPoints), ...
    'snapshotPoints', numel(snapshotPoints), ...
    'methods', methodCount, ...
    'fixedSnapshots', paper.fixedSnapshotCount, ...
    'fixedSnrDb', paper.fixedSnrDb);
result = uav.experiment.finalize(result, summaryTable, fig, cfg, ...
    actualCounts, timerValue);
end

function summary = summarizeSweep(values, cfg, experimentIndex, pointOffset)
%SUMMARIZESWEEP 将[试次×扫描点×方法]逐点汇总为[扫描点×方法]。
pointCount = size(values, 2);
methodCount = size(values, 3);
fields = {'medianLinear', 'lowLinear', 'meanLinear', 'medianDb', ...
    'lowDb', 'meanDb'};
for fieldIndex = 1:numel(fields)
    summary.(fields{fieldIndex}) = nan(pointCount, methodCount);
end
summary.ci95Linear = nan(pointCount, 2, methodCount);
summary.ci95Db = nan(pointCount, 2, methodCount);
summary.validCounts = zeros(pointCount, methodCount);
for pointIndex = 1:pointCount
    samples = reshape(values(:, pointIndex, :), size(values, 1), methodCount);
    bootstrapSubstream = uav.experiment.substream(cfg, experimentIndex, ...
        pointOffset + pointIndex, 0, 0, "bootstrap");
    point = uav.experiment.summarizeLinear(samples, cfg, bootstrapSubstream);
    for fieldIndex = 1:numel(fields)
        summary.(fields{fieldIndex})(pointIndex, :) = point.(fields{fieldIndex});
    end
    summary.ci95Linear(pointIndex, :, :) = reshape(point.ci95Linear, 1, 2, []);
    summary.ci95Db(pointIndex, :, :) = reshape(point.ci95Db, 1, 2, []);
    summary.validCounts(pointIndex, :) = point.validCounts;
end
summary.aggregationDomain = "linear";
summary.primaryStatistic = "mean";
end

function fig = plotReproduction(snrPoints, snapshotPoints, methods, summary, cfg)
%PLOTREPRODUCTION 绘制与论文图1(a)、图1(c)相对应的两幅曲线。
fig = uav.plot.newFigure(cfg);
% 两个子图共用用户设置的25号字体时，按图高派生最小宽度，避免标题和图例重叠。
position = fig.Position;
position(3) = max(position(3), 2 * position(4));
fig.Position = position;
labels = methodLabels(methods);
styles = ["--o", "-s", "-.^", "-d"];
colors = lines(numel(methods));

subplot(1, 2, 1);
hold on;
for methodIndex = 1:numel(methods)
    plot(snrPoints, summary.snr.meanDb(:, methodIndex), ...
        styles(methodIndex), 'LineWidth', cfg.plot.lineWidth, ...
        'MarkerSize', cfg.plot.markerSize, 'Color', colors(methodIndex, :));
end
xlabel('输入SNR/dB', 'FontSize', cfg.plot.fontSize);
ylabel('平均输出SINR/dB', 'FontSize', cfg.plot.fontSize);
title('Gu–Leshem 2012例1：SNR扫描', 'FontSize', cfg.plot.fontSize);
legend(labels, 'Location', 'best', 'FontSize', cfg.plot.fontSize);
set(gca, 'FontSize', cfg.plot.fontSize); grid on; box on;

subplot(1, 2, 2);
hold on;
for methodIndex = 1:numel(methods)
    plot(snapshotPoints, summary.snapshot.meanDb(:, methodIndex), ...
        styles(methodIndex), 'LineWidth', cfg.plot.lineWidth, ...
        'MarkerSize', cfg.plot.markerSize, 'Color', colors(methodIndex, :));
end
xlabel('快拍数K', 'FontSize', cfg.plot.fontSize);
ylabel('平均输出SINR/dB', 'FontSize', cfg.plot.fontSize);
title('Gu–Leshem 2012例1：快拍数扫描', 'FontSize', cfg.plot.fontSize);
legend(labels, 'Location', 'best', 'FontSize', cfg.plot.fontSize);
set(gca, 'FontSize', cfg.plot.fontSize); grid on; box on;
end

function labels = methodLabels(methods)
known = ["oracle", "esp", "wco", "gu_leshem2012_rec"];
names = ["理想先验MVDR", "特征空间法", "最坏情况法", "Gu–Leshem Rec"];
[found, indices] = ismember(methods, known);
if ~all(found)
    error('uav:experiment:UnknownGuLeshemLabel', ...
        'Every paper comparison method must have a Chinese label.');
end
labels = names(indices);
end

function output = reproductionTable(snrPoints, snapshotPoints, methods, summary)
%REPRODUCTIONTABLE 生成两类扫描共用的长表CSV。
rowCount = (numel(snrPoints) + numel(snapshotPoints)) * numel(methods);
sweepType = strings(rowCount, 1);
xValue = zeros(rowCount, 1);
method = strings(rowCount, 1);
meanSinrDb = zeros(rowCount, 1);
medianSinrDb = zeros(rowCount, 1);
row = 0;
for methodIndex = 1:numel(methods)
    for pointIndex = 1:numel(snrPoints)
        row = row + 1;
        sweepType(row) = "snr_db";
        xValue(row) = snrPoints(pointIndex);
        method(row) = methods(methodIndex);
        meanSinrDb(row) = summary.snr.meanDb(pointIndex, methodIndex);
        medianSinrDb(row) = summary.snr.medianDb(pointIndex, methodIndex);
    end
    for pointIndex = 1:numel(snapshotPoints)
        row = row + 1;
        sweepType(row) = "snapshot_count";
        xValue(row) = snapshotPoints(pointIndex);
        method(row) = methods(methodIndex);
        meanSinrDb(row) = summary.snapshot.meanDb(pointIndex, methodIndex);
        medianSinrDb(row) = summary.snapshot.medianDb(pointIndex, methodIndex);
    end
end
output = table(sweepType, xValue, method, meanSinrDb, medianSinrDb);
end

function validateConfiguration(cfg)
if ~isstruct(cfg) || ~isfield(cfg, 'literature') || ...
        ~isfield(cfg.literature, 'guLeshem2012')
    error('uav:experiment:MissingGuLeshemConfiguration', ...
        'cfg.literature.guLeshem2012 is required.');
end
paper = cfg.literature.guLeshem2012;
if paper.fixedSnapshotCount < paper.sensorCount
    error('uav:experiment:SingularPaperScenario', ...
        'The fixed paper reproduction requires at least M snapshots.');
end
if any(paper.snapshotSweep < paper.sensorCount)
    error('uav:experiment:SingularPaperScenario', ...
        'Snapshot-sweep values must be at least M for the unmodified Rec formula.');
end
end
