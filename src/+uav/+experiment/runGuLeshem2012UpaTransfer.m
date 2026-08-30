function result = runGuLeshem2012UpaTransfer(cfg)
%RUNGULESHEM2012UPATRANSFER 评估Rec二维直接扩展的网格和扇区敏感性。
%   输出SINR、INCM NMSE、主特征值比例、干扰子空间距离，以及论文核心
%   扇区假设在姿态抖动下的违约率。真值只用于Oracle和实验后评价。

timerValue = tic;
% 入口重新校验，确保方法集合、快拍数和二维网格满足文献专项合同。
cfg = uav.config.validate(cfg);
paper = cfg.literature.guLeshem2012;
methods = reshape(string(paper.methodNames), 1, []);
jitterPoints = reshape(paper.upaTransferJitterDeg, 1, []);
gridSizes = reshape(paper.upaTransferGridSizes, 1, []);
trialCount = paper.monteCarloTrials;
periodCount = paper.upaTransferPeriods;
methodCount = numel(methods);
J = cfg.derived.J;
validateTransfer(cfg, jitterPoints, gridSizes, periodCount);

raw = struct( ...
    'sinrLinear', zeros(trialCount, periodCount, numel(jitterPoints), ...
        numel(gridSizes), methodCount), ...
    'incmNmse', nan(trialCount, periodCount, numel(jitterPoints), numel(gridSizes)), ...
    'largestEigenvalueRatio', nan(trialCount, periodCount, ...
        numel(jitterPoints), numel(gridSizes)), ...
    'interferenceSubspaceDistance', nan(trialCount, periodCount, ...
        numel(jitterPoints), numel(gridSizes)), ...
    'normalizedEigenvalueProfileError', nan(trialCount, periodCount, ...
        numel(jitterPoints), numel(gridSizes)), ...
    'estimatedConditionNumber', nan(trialCount, periodCount, ...
        numel(jitterPoints), numel(gridSizes)), ...
    'desiredOutsideExclusion', false(trialCount, periodCount, ...
        numel(jitterPoints), numel(gridSizes)), ...
    'interfererInsideExclusion', false(trialCount, periodCount, ...
        numel(jitterPoints), numel(gridSizes), J));
quadratureDiagnostics = cell(numel(gridSizes), 1);

for gridIndex = 1:numel(gridSizes)
    for jitterIndex = 1:numel(jitterPoints)
        pointCfg = cfg;
        pointCfg.spectrum.gridSize = gridSizes(gridIndex);
        pointCfg.attitude.stdDeg = repmat(jitterPoints(jitterIndex), 1, 3);
        pointCfg.timeline.L = paper.upaTransferSnapshotCount;
        pointCfg = uav.config.validate(pointCfg);
        [scanManifold, scanWeights, quadratureDiagnostics{gridIndex}] = ...
            uav.spectrum.guLeshem2012UpaQuadrature(pointCfg);
        for trialIndex = 1:trialCount
            % 网格大小不进入随机编号，保证不同网格使用逐元素相同的场景和快拍。
            scenarioSubstream = uav.experiment.substream(cfg, 8, ...
                jitterIndex, trialIndex, 0, "scenario");
            scenario = uav.scenario.generate(pointCfg, ...
                uav.util.makeStream(cfg.random.baseSeed, scenarioSubstream));
            for periodIndex = 1:periodCount
                snapshotSubstream = uav.experiment.substream(cfg, 8, ...
                    jitterIndex, trialIndex, periodIndex, "snapshot");
                obs = uav.scenario.generatePeriod(pointCfg, scenario, periodIndex, ...
                    uav.util.makeStream(cfg.random.baseSeed, snapshotSubstream));
                safe = uav.experiment.algorithmObservation(obs);
                adapted = struct( ...
                    'X', safe.X, ...
                    'Rx', safe.Rx, ...
                    'nominalSteering', pointCfg.derived.nominalSteeringReference, ...
                    'scanManifold', scanManifold, ...
                    'scanWeights', scanWeights, ...
                    'truth', obs.truth);
                record = uav.experiment.runGuLeshem2012Methods( ...
                    adapted, pointCfg, methods);
                raw.sinrLinear(trialIndex, periodIndex, jitterIndex, ...
                    gridIndex, :) = reshape(record.sinrLinear, 1, 1, 1, 1, []);

                recIndex = find(methods == "gu_leshem2012_rec", 1);
                estimatedIncm = record.diagnostics{recIndex}.estimatedIncm;
                raw.incmNmse(trialIndex, periodIndex, jitterIndex, gridIndex) = ...
                    record.incmNmse(recIndex);
                raw.largestEigenvalueRatio(trialIndex, periodIndex, ...
                    jitterIndex, gridIndex) = largestEigenvalueRatio( ...
                    estimatedIncm, obs.truth.trueIncm);
                raw.interferenceSubspaceDistance(trialIndex, periodIndex, ...
                    jitterIndex, gridIndex) = subspaceDistance( ...
                    estimatedIncm, obs.truth.trueIncm, J);
                raw.normalizedEigenvalueProfileError(trialIndex, periodIndex, ...
                    jitterIndex, gridIndex) = normalizedEigenvalueProfileError( ...
                    estimatedIncm, obs.truth.trueIncm);
                raw.estimatedConditionNumber(trialIndex, periodIndex, ...
                    jitterIndex, gridIndex) = cond(estimatedIncm);
                [desiredOutside, interfererInside] = sectorViolations(obs, pointCfg);
                raw.desiredOutsideExclusion(trialIndex, periodIndex, ...
                    jitterIndex, gridIndex) = desiredOutside;
                raw.interfererInsideExclusion(trialIndex, periodIndex, ...
                    jitterIndex, gridIndex, :) = reshape(interfererInside, 1, 1, 1, 1, []);
            end
        end
    end
end

summary = summarizeTransfer(raw, cfg, jitterPoints, gridSizes, methods);
fig = plotTransfer(jitterPoints, gridSizes, methods, summary, cfg);
summaryTable = transferTable(jitterPoints, gridSizes, methods, summary);
result = struct( ...
    'experimentName', "exp08_gu_leshem2012_upa_transfer", ...
    'methodNames', methods, ...
    'jitterPointsDeg', jitterPoints, ...
    'gridSizes', gridSizes, ...
    'quadratureDiagnostics', {quadratureDiagnostics}, ...
    'raw', raw, ...
    'summary', summary);
actualCounts = struct( ...
    'trials', trialCount, 'periods', periodCount, ...
    'jitterPoints', numel(jitterPoints), ...
    'gridSizes', numel(gridSizes), 'methods', methodCount, ...
    'snapshots', paper.upaTransferSnapshotCount);
result = uav.experiment.finalize(result, summaryTable, fig, cfg, ...
    actualCounts, timerValue);
end

function summary = summarizeTransfer(raw, cfg, jitterPoints, gridSizes, methods)
%SUMMARIZETRANSFER 以试次-周期样本为行，逐抖动点和网格汇总。
summary.sinr.meanDb = nan(numel(jitterPoints), numel(gridSizes), numel(methods));
summary.sinr.medianDb = nan(size(summary.sinr.meanDb));
summary.incmNmseMean = nan(numel(jitterPoints), numel(gridSizes));
summary.largestEigenvalueRatioMean = nan(numel(jitterPoints), numel(gridSizes));
summary.interferenceSubspaceDistanceMean = nan(numel(jitterPoints), numel(gridSizes));
summary.normalizedEigenvalueProfileErrorMean = ...
    nan(numel(jitterPoints), numel(gridSizes));
summary.estimatedConditionNumberMedian = ...
    nan(numel(jitterPoints), numel(gridSizes));
summary.desiredOutsideRate = nan(numel(jitterPoints), numel(gridSizes));
summary.interfererInsideRate = nan(numel(jitterPoints), numel(gridSizes), ...
    size(raw.interfererInsideExclusion, 5));
for jitterIndex = 1:numel(jitterPoints)
    for gridIndex = 1:numel(gridSizes)
        samples = reshape(raw.sinrLinear(:, :, jitterIndex, gridIndex, :), ...
            [], numel(methods));
        bootstrapSubstream = uav.experiment.substream(cfg, 8, ...
            (gridIndex - 1) * numel(jitterPoints) + jitterIndex, ...
            0, 0, "bootstrap");
        point = uav.experiment.summarizeLinear(samples, cfg, bootstrapSubstream);
        summary.sinr.meanDb(jitterIndex, gridIndex, :) = ...
            reshape(point.meanDb, 1, 1, []);
        summary.sinr.medianDb(jitterIndex, gridIndex, :) = ...
            reshape(point.medianDb, 1, 1, []);
        summary.incmNmseMean(jitterIndex, gridIndex) = mean( ...
            raw.incmNmse(:, :, jitterIndex, gridIndex), 'all');
        summary.largestEigenvalueRatioMean(jitterIndex, gridIndex) = mean( ...
            raw.largestEigenvalueRatio(:, :, jitterIndex, gridIndex), 'all');
        summary.interferenceSubspaceDistanceMean(jitterIndex, gridIndex) = mean( ...
            raw.interferenceSubspaceDistance(:, :, jitterIndex, gridIndex), 'all');
        summary.normalizedEigenvalueProfileErrorMean(jitterIndex, gridIndex) = mean( ...
            raw.normalizedEigenvalueProfileError(:, :, jitterIndex, gridIndex), 'all');
        summary.estimatedConditionNumberMedian(jitterIndex, gridIndex) = median( ...
            raw.estimatedConditionNumber(:, :, jitterIndex, gridIndex), 'all');
        summary.desiredOutsideRate(jitterIndex, gridIndex) = mean( ...
            raw.desiredOutsideExclusion(:, :, jitterIndex, gridIndex), 'all');
        for interfererIndex = 1:size(raw.interfererInsideExclusion, 5)
            summary.interfererInsideRate(jitterIndex, gridIndex, interfererIndex) = ...
                mean(raw.interfererInsideExclusion(:, :, jitterIndex, ...
                gridIndex, interfererIndex), 'all');
        end
    end
end
end

function ratio = largestEigenvalueRatio(estimate, truth)
ratio = max(real(eig((estimate + estimate') / 2))) / ...
    max(real(eig((truth + truth') / 2)));
end

function distance = subspaceDistance(estimate, truth, dimension)
if dimension == 0
    distance = 0;
    return;
end

[estimateVectors, estimateValues] = eig((estimate + estimate') / 2, 'vector');
[~, estimateOrder] = sort(real(estimateValues), 'descend');
[truthVectors, truthValues] = eig((truth + truth') / 2, 'vector');
[~, truthOrder] = sort(real(truthValues), 'descend');
estimateBasis = estimateVectors(:, estimateOrder(1:dimension));
truthBasis = truthVectors(:, truthOrder(1:dimension));
distance = norm(estimateBasis * estimateBasis' - ...
    truthBasis * truthBasis', 'fro');
end

function errorValue = normalizedEigenvalueProfileError(estimate, truth)
%NORMALIZEDEIGENVALUEPROFILEERROR 比较去除整体尺度后的特征值能量分布。
%   MVDR对协方差的统一正比例缩放不敏感，因此该指标用于把“绝对功率
%   尺度错误”与“相对特征值结构错误”分开，避免只凭最大特征值下结论。
estimateValues = sort(max(real(eig((estimate + estimate') / 2)), 0), 'descend');
truthValues = sort(max(real(eig((truth + truth') / 2)), 0), 'descend');
estimateTotal = sum(estimateValues);
truthTotal = sum(truthValues);
if estimateTotal <= 0 || truthTotal <= 0
    error('uav:experiment:InvalidIncmEigenvalueProfile', ...
        'INCM eigenvalue profiles must have positive total power.');
end
errorValue = norm(estimateValues / estimateTotal - ...
    truthValues / truthTotal, 2);
end

function [desiredOutside, interfererInside] = sectorViolations(obs, cfg)
nominalDesired = cfg.derived.nominalArrayUnitVectors(:, 1);
actualDirections = obs.truth.arrayUnitVectors;
cosines = max(-1, min(1, real(nominalDesired' * actualDirections)));
distances = acos(cosines);
halfWidth = deg2rad(cfg.spectrum.desiredExclusionHalfWidthDeg);
desiredOutside = distances(1) > halfWidth;
interfererInside = distances(2:end) <= halfWidth;
end

function fig = plotTransfer(jitterPoints, gridSizes, methods, summary, cfg)
fig = uav.plot.newFigure(cfg);
labels = methodLabels(methods);
styles = ["--o", "-s", "-.^", "-d"];
colors = lines(numel(methods));
panelCount = numel(gridSizes) + 1;
% 面板数随网格配置变化，按图高自动派生总宽度，避免大字号标题相互覆盖。
position = fig.Position;
position(3) = max(position(3), panelCount * position(4));
fig.Position = position;
for gridIndex = 1:numel(gridSizes)
    subplot(1, panelCount, gridIndex); hold on;
    for methodIndex = 1:numel(methods)
        plot(jitterPoints, summary.sinr.meanDb(:, gridIndex, methodIndex), ...
            styles(methodIndex), 'LineWidth', cfg.plot.lineWidth, ...
            'MarkerSize', cfg.plot.markerSize, 'Color', colors(methodIndex, :));
    end
    xlabel('三轴姿态抖动RMS/°', 'FontSize', cfg.plot.fontSize);
    ylabel('平均输出SINR/dB', 'FontSize', cfg.plot.fontSize);
    title(sprintf('%d×%d方向余弦网格', gridSizes(gridIndex), ...
        gridSizes(gridIndex)), 'FontSize', cfg.plot.fontSize);
    legend(labels, 'Location', 'best', 'FontSize', cfg.plot.fontSize);
    set(gca, 'FontSize', cfg.plot.fontSize); grid on; box on;
end
subplot(1, panelCount, panelCount); hold on;
for gridIndex = 1:numel(gridSizes)
    plot(jitterPoints, summary.incmNmseMean(:, gridIndex), ...
        '-o', 'LineWidth', cfg.plot.lineWidth, ...
        'MarkerSize', cfg.plot.markerSize);
end
xlabel('三轴姿态抖动RMS/°', 'FontSize', cfg.plot.fontSize);
ylabel('INCM归一化平方误差', 'FontSize', cfg.plot.fontSize);
title('二维直接扩展的网格敏感性', 'FontSize', cfg.plot.fontSize);
gridLabels = string(gridSizes) + "×" + string(gridSizes);
legend(gridLabels, ...
    'Location', 'best', 'FontSize', cfg.plot.fontSize);
set(gca, 'FontSize', cfg.plot.fontSize); grid on; box on;
end

function labels = methodLabels(methods)
%METHODLABELS 按用户配置的方法顺序返回对应中文图例。
known = ["oracle", "esp", "wco", "gu_leshem2012_rec"];
names = ["理想先验MVDR", "特征空间法", "最坏情况法", "Gu–Leshem Rec"];
[found, indices] = ismember(methods, known);
if ~all(found)
    error('uav:experiment:UnknownGuLeshemLabel', ...
        'Every paper comparison method must have a Chinese label.');
end
labels = names(indices);
end

function output = transferTable(jitterPoints, gridSizes, methods, summary)
rowCount = numel(jitterPoints) * numel(gridSizes) * numel(methods);
jitterRmsDeg = zeros(rowCount, 1);
gridSize = zeros(rowCount, 1);
method = strings(rowCount, 1);
meanSinrDb = zeros(rowCount, 1);
recIncmNmseMean = zeros(rowCount, 1);
sharedDesiredOutsideRate = zeros(rowCount, 1);
sharedMaximumInterfererInsideRate = zeros(rowCount, 1);
recNormalizedEigenvalueProfileErrorMean = zeros(rowCount, 1);
recEstimatedConditionNumberMedian = zeros(rowCount, 1);
row = 0;
for gridIndex = 1:numel(gridSizes)
    for jitterIndex = 1:numel(jitterPoints)
        for methodIndex = 1:numel(methods)
            row = row + 1;
            jitterRmsDeg(row) = jitterPoints(jitterIndex);
            gridSize(row) = gridSizes(gridIndex);
            method(row) = methods(methodIndex);
            meanSinrDb(row) = summary.sinr.meanDb( ...
                jitterIndex, gridIndex, methodIndex);
            recIncmNmseMean(row) = summary.incmNmseMean(jitterIndex, gridIndex);
            sharedDesiredOutsideRate(row) = ...
                summary.desiredOutsideRate(jitterIndex, gridIndex);
            sharedMaximumInterfererInsideRate(row) = max( ...
                summary.interfererInsideRate(jitterIndex, gridIndex, :), [], 'all');
            recNormalizedEigenvalueProfileErrorMean(row) = ...
                summary.normalizedEigenvalueProfileErrorMean(jitterIndex, gridIndex);
            recEstimatedConditionNumberMedian(row) = ...
                summary.estimatedConditionNumberMedian(jitterIndex, gridIndex);
        end
    end
end
output = table(jitterRmsDeg, gridSize, method, meanSinrDb, ...
    recIncmNmseMean, recNormalizedEigenvalueProfileErrorMean, ...
    recEstimatedConditionNumberMedian, sharedDesiredOutsideRate, ...
    sharedMaximumInterfererInsideRate);
end

function validateTransfer(cfg, jitterPoints, gridSizes, periodCount)
if cfg.literature.guLeshem2012.upaTransferSnapshotCount < cfg.derived.M
    error('uav:experiment:SingularUpaTransferScenario', ...
        'UPA transfer requires at least M snapshots for the unmodified Rec formula.');
end
if isempty(jitterPoints) || any(~isfinite(jitterPoints)) || any(jitterPoints < 0)
    error('uav:experiment:InvalidUpaTransferJitter', ...
        'UPA transfer jitter points must be finite and nonnegative.');
end
if isempty(gridSizes) || any(gridSizes < 3) || any(mod(gridSizes, 2) ~= 1)
    error('uav:experiment:InvalidUpaTransferGrid', ...
        'UPA transfer grid sizes must be odd integers at least three.');
end
if periodCount < 1 || periodCount > cfg.timeline.K || periodCount ~= floor(periodCount)
    error('uav:experiment:InvalidUpaTransferPeriods', ...
        'UPA transfer periods must be a positive integer no larger than K.');
end
end
