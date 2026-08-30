function result = runUpaIncmDiagnosis(cfg)
%RUNUPAINCMDIAGNOSIS 编排UPA二维INCM诊断实验9并导出完整结果包。
%   零抖动分支复用每个试次唯一的64快拍观测，8快拍严格取其前缀；
%   抖动分支按试次、抖动点和周期配对两种快拍。真实量只在权值和INCM
%   已形成后用于指标评价，绝不进入候选、积分、拟合或权值参数选择。

timerValue = tic;
cfg = uav.config.validate(cfg);
diagnosis = cfg.literature.guLeshem2012.upaDiagnosis;
methodNames = ["oracle", "esp", "wco", "gu_leshem2012_rec", ...
    "upa_covfit_incm"];
snapshotCounts = reshape(diagnosis.snapshotCounts, 1, []);
gridSizes = reshape(diagnosis.gridSizes, 1, []);
measureNames = reshape(string(diagnosis.measureNames), 1, []);
jitterPoints = reshape(diagnosis.jitterSweepDeg, 1, []);

zeroJitter = collectZeroJitter(cfg, diagnosis, methodNames, ...
    snapshotCounts, gridSizes, measureNames);
jitterSweep = collectJitterSweep(cfg, diagnosis, methodNames, ...
    snapshotCounts, jitterPoints);
summary = summarizeDiagnosis(zeroJitter, jitterSweep, cfg, diagnosis, ...
    snapshotCounts, gridSizes, measureNames, jitterPoints, methodNames);
fig = plotDiagnosis(summary, cfg, snapshotCounts, gridSizes, ...
    measureNames, jitterPoints, methodNames, diagnosis);
summaryTable = diagnosisTable(summary, snapshotCounts, gridSizes, ...
    measureNames, jitterPoints, methodNames);

result = struct( ...
    'experimentName', "exp09_upa_incm_diagnosis", ...
    'methodNames', methodNames, ...
    'snapshotCounts', snapshotCounts, ...
    'gridSizes', gridSizes, ...
    'measureNames', measureNames, ...
    'jitterPointsDeg', jitterPoints, ...
    'raw', struct('zeroJitter', zeroJitter, 'jitterSweep', jitterSweep), ...
    'summary', summary);
actualCounts = struct( ...
    'trials', diagnosis.monteCarloTrials, ...
    'periods', diagnosis.periods, ...
    'zeroJitterPeriods', 1, ...
    'snapshotCounts', snapshotCounts, ...
    'gridSizes', gridSizes, ...
    'measureNames', measureNames, ...
    'jitterPoints', numel(jitterPoints), ...
    'methods', numel(methodNames));
result = uav.experiment.finalize(result, summaryTable, fig, cfg, ...
    actualCounts, timerValue);
end

function raw = collectZeroJitter(cfg, diagnosis, methods, snapshots, grids, measures)
%COLLECTZEROJITTER 对同一零抖动快拍执行全部网格、测度和方法。
trialCount = diagnosis.monteCarloTrials;
methodCount = numel(methods);
J = cfg.derived.J;
raw = struct( ...
    'sinrLinear', nan(trialCount, numel(snapshots), numel(grids), ...
        numel(measures), methodCount), ...
    'status', repmat("", trialCount, numel(snapshots), numel(grids), ...
        numel(measures), methodCount), ...
    'reasonCode', repmat("", trialCount, numel(snapshots), numel(grids), ...
        numel(measures), methodCount), ...
    'snapshotPrefixMatched', false(trialCount, numel(grids), numel(measures)), ...
    'desiredOutsideExclusion', false(trialCount, numel(grids), numel(measures)), ...
    'interfererInsideExclusion', false(trialCount, numel(grids), ...
        numel(measures), J), ...
    'originalRec', emptyIncmDiagnostics( ...
        [trialCount, numel(snapshots), numel(grids), numel(measures)], J), ...
    'upaCovfit', emptyIncmDiagnostics( ...
        [trialCount, numel(snapshots), numel(grids), numel(measures)], J));

generationCfg = cfg;
generationCfg.timeline.L = max(snapshots);
generationCfg.attitude.stdDeg = [0, 0, 0];
generationCfg = uav.config.validate(generationCfg);
for trialIndex = 1:trialCount
    scenarioStream = uav.experiment.substream(cfg, diagnosis.experimentIndex, ...
        1, trialIndex, 0, "scenario");
    scenario = uav.scenario.generate(generationCfg, ...
        uav.util.makeStream(cfg.random.baseSeed, scenarioStream));
    snapshotStream = uav.experiment.substream(cfg, diagnosis.experimentIndex, ...
        1, trialIndex, 1, "snapshot");
    obs64 = uav.scenario.generatePeriod(generationCfg, scenario, 1, ...
        uav.util.makeStream(cfg.random.baseSeed, snapshotStream));
    paired = cell(1, numel(snapshots));
    for snapshotIndex = 1:numel(snapshots)
        paired{snapshotIndex} = uav.scenario.truncatePeriodObservation( ...
            obs64, snapshots(snapshotIndex));
    end
    prefixMatched = isequal(paired{1}.X, ...
        paired{end}.X(:, 1:snapshots(1)));

    for gridIndex = 1:numel(grids)
        for measureIndex = 1:numel(measures)
            for snapshotIndex = 1:numel(snapshots)
                % 每个点重新校验，以按当前网格刷新加载谱使用的流形缓存。
                pointCfg = generationCfg;
                pointCfg.timeline.L = snapshots(snapshotIndex);
                pointCfg.spectrum.gridSize = grids(gridIndex);
                pointCfg = uav.config.validate(pointCfg);
                integrationGrid = uav.spectrum.upaIntegrationGrid( ...
                    pointCfg, grids(gridIndex), measures(measureIndex));
                obs = paired{snapshotIndex};
                record = uav.experiment.runUpaDiagnosisMethods( ...
                    obs, pointCfg, methods, integrationGrid);
                raw.sinrLinear(trialIndex, snapshotIndex, gridIndex, ...
                    measureIndex, :) = reshape(record.sinrLinear, 1, 1, 1, 1, []);
                raw.status(trialIndex, snapshotIndex, gridIndex, ...
                    measureIndex, :) = reshape(record.status, 1, 1, 1, 1, []);
                raw.reasonCode(trialIndex, snapshotIndex, gridIndex, ...
                    measureIndex, :) = reshape(record.reasonCode, 1, 1, 1, 1, []);
                raw.snapshotPrefixMatched(trialIndex, gridIndex, measureIndex) = ...
                    prefixMatched;
                [desiredOutside, interfererInside] = sectorViolations(obs, pointCfg);
                raw.desiredOutsideExclusion(trialIndex, gridIndex, measureIndex) = ...
                    desiredOutside;
                if J > 0
                    raw.interfererInsideExclusion(trialIndex, gridIndex, ...
                        measureIndex, :) = reshape(interfererInside, 1, 1, 1, []);
                end
                raw = storeMethodDiagnostics(raw, record, methods, obs, ...
                    integrationGrid, pointCfg, ...
                    [trialIndex, snapshotIndex, gridIndex, measureIndex]);
            end
        end
    end
end
end

function raw = collectJitterSweep(cfg, diagnosis, methods, snapshots, jitterPoints)
%COLLECTJITTERSWEEP 对选定网格和测度执行动态配对扫描。
trialCount = diagnosis.monteCarloTrials;
periodCount = diagnosis.periods;
methodCount = numel(methods);
J = cfg.derived.J;
baseDimensions = [trialCount, periodCount, numel(jitterPoints), numel(snapshots)];
raw = struct( ...
    'sinrLinear', nan(trialCount, periodCount, numel(jitterPoints), ...
        numel(snapshots), methodCount), ...
    'status', repmat("", trialCount, periodCount, numel(jitterPoints), ...
        numel(snapshots), methodCount), ...
    'reasonCode', repmat("", trialCount, periodCount, numel(jitterPoints), ...
        numel(snapshots), methodCount), ...
    'snapshotPrefixMatched', false(trialCount, periodCount, numel(jitterPoints)), ...
    'upaCovfitWcoSinrLinear', nan(trialCount, periodCount, ...
        numel(jitterPoints), numel(snapshots)), ...
    'upaCovfitWcoStatus', repmat("", trialCount, periodCount, ...
        numel(jitterPoints), numel(snapshots)), ...
    'desiredOutsideExclusion', false(trialCount, periodCount, ...
        numel(jitterPoints)), ...
    'interfererInsideExclusion', false(trialCount, periodCount, ...
        numel(jitterPoints), J), ...
    'originalRec', emptyIncmDiagnostics(baseDimensions, J), ...
    'upaCovfit', emptyIncmDiagnostics(baseDimensions, J));

for jitterIndex = 1:numel(jitterPoints)
    generationCfg = cfg;
    generationCfg.timeline.L = max(snapshots);
    generationCfg.spectrum.gridSize = diagnosis.selectedGridSize;
    generationCfg.attitude.stdDeg = repmat(jitterPoints(jitterIndex), 1, 3);
    generationCfg = uav.config.validate(generationCfg);
    for trialIndex = 1:trialCount
        scenarioStream = uav.experiment.substream(cfg, diagnosis.experimentIndex, ...
            jitterIndex + 1, trialIndex, 0, "scenario");
        scenario = uav.scenario.generate(generationCfg, ...
            uav.util.makeStream(cfg.random.baseSeed, scenarioStream));
        for periodIndex = 1:periodCount
            snapshotStream = uav.experiment.substream(cfg, ...
                diagnosis.experimentIndex, jitterIndex + 1, trialIndex, ...
                periodIndex, "snapshot");
            obs64 = uav.scenario.generatePeriod(generationCfg, scenario, ...
                periodIndex, uav.util.makeStream(cfg.random.baseSeed, snapshotStream));
            paired = cell(1, numel(snapshots));
            for snapshotIndex = 1:numel(snapshots)
                paired{snapshotIndex} = uav.scenario.truncatePeriodObservation( ...
                    obs64, snapshots(snapshotIndex));
            end
            raw.snapshotPrefixMatched(trialIndex, periodIndex, jitterIndex) = ...
                isequal(paired{1}.X, paired{end}.X(:, 1:snapshots(1)));

            for snapshotIndex = 1:numel(snapshots)
                % 快拍数变化也重新校验，避免沿用L=64的派生谱缓存。
                pointCfg = generationCfg;
                pointCfg.timeline.L = snapshots(snapshotIndex);
                pointCfg = uav.config.validate(pointCfg);
                integrationGrid = uav.spectrum.upaIntegrationGrid(pointCfg, ...
                    diagnosis.selectedGridSize, diagnosis.selectedMeasureName);
                obs = paired{snapshotIndex};
                record = uav.experiment.runUpaDiagnosisMethods( ...
                    obs, pointCfg, methods, integrationGrid);
                raw.sinrLinear(trialIndex, periodIndex, jitterIndex, ...
                    snapshotIndex, :) = reshape(record.sinrLinear, 1, 1, 1, 1, []);
                raw.status(trialIndex, periodIndex, jitterIndex, ...
                    snapshotIndex, :) = reshape(record.status, 1, 1, 1, 1, []);
                raw.reasonCode(trialIndex, periodIndex, jitterIndex, ...
                    snapshotIndex, :) = reshape(record.reasonCode, 1, 1, 1, 1, []);

                upaIndex = find(methods == "upa_covfit_incm", 1);
                estimatedIncm = record.diagnostics{upaIndex}.estimatedIncm;
                methodCfg = pointCfg;
                methodCfg.beamforming.epsilonS = sqrt( ...
                    cfg.literature.guLeshem2012.wcoSquaredRadiusScale * ...
                    cfg.derived.M);
                [wcoWeight, ~] = uav.beamforming.solve("wco", estimatedIncm, ...
                    obs.nominalSteeringReference, methodCfg);
                wcoMetrics = uav.metrics.evaluate(wcoWeight, obs.truth, methodCfg);
                raw.upaCovfitWcoSinrLinear(trialIndex, periodIndex, ...
                    jitterIndex, snapshotIndex) = wcoMetrics.sinrLinear;
                raw.upaCovfitWcoStatus(trialIndex, periodIndex, ...
                    jitterIndex, snapshotIndex) = "ok";
                % 五方法与独立WCO权值均已形成后，才允许读取真值诊断。
                raw = storeMethodDiagnostics(raw, record, methods, obs, ...
                    integrationGrid, pointCfg, ...
                    [trialIndex, periodIndex, jitterIndex, snapshotIndex]);
            end
            % 扇区违约与INCM误差一样仅是后验评价，不参与下一周期计算。
            [desiredOutside, interfererInside] = sectorViolations(obs64, generationCfg);
            raw.desiredOutsideExclusion(trialIndex, periodIndex, jitterIndex) = ...
                desiredOutside;
            if J > 0
                raw.interfererInsideExclusion(trialIndex, periodIndex, ...
                    jitterIndex, :) = reshape(interfererInside, 1, 1, 1, []);
            end
        end
    end
end
end

function raw = storeMethodDiagnostics(raw, record, methods, obs, grid, pointCfg, indices)
%STOREMETHODDIAGNOSTICS 在算法结束后才读取真值并写入两类独立诊断数组。
recIndex = find(methods == "gu_leshem2012_rec", 1);
upaIndex = find(methods == "upa_covfit_incm", 1);
gridSize = pointCfg.spectrum.gridSize;
radiusBins = pointCfg.literature.guLeshem2012.upaDiagnosis.peakNmsRadiusBins;
if record.status(recIndex) == "ok"
    recDiagnostics = record.diagnostics{recIndex};
    metrics = incmMetrics(recDiagnostics.estimatedIncm, ...
        obs.truth.trueIncm, size(obs.truth.trueASVs, 2) - 1);
    local = localRecDiagnostics(recDiagnostics.reconstruction, grid, ...
        obs.truth.directionCosines(:, 2:end).', gridSize, radiusBins);
    raw.originalRec = assignDiagnostics(raw.originalRec, indices, metrics, local);
end
upaDiagnostics = record.diagnostics{upaIndex};
metrics = incmMetrics(upaDiagnostics.estimatedIncm, obs.truth.trueIncm, ...
    size(obs.truth.trueASVs, 2) - 1);
local = localCovfitDiagnostics(upaDiagnostics.reconstruction, ...
    obs.truth.directionCosines(:, 2:end).', gridSize, radiusBins);
raw.upaCovfit = assignDiagnostics(raw.upaCovfit, indices, metrics, local);
end

function diagnostics = emptyIncmDiagnostics(baseDimensions, interfererCount)
%EMPTYINCMDIAGNOSTICS 以NaN统一表示方法不适用或局部真值邻域无候选。
diagnostics = struct( ...
    'incmNmse', nan(baseDimensions), ...
    'traceRatio', nan(baseDimensions), ...
    'largestEigenvalueRatio', nan(baseDimensions), ...
    'normalizedEigenvalueProfileError', nan(baseDimensions), ...
    'conditionNumber', nan(baseDimensions), ...
    'interferenceSubspaceDistance', nan(baseDimensions), ...
    'localPeakValue', nan([baseDimensions, interfererCount]), ...
    'localIntegratedEnergy', nan([baseDimensions, interfererCount]), ...
    'localPeakWidth', nan([baseDimensions, interfererCount]));
end

function output = assignDiagnostics(output, indices, metrics, local)
%ASSIGNDIAGNOSTICS 按统一索引写入标量和逐干扰后验指标。
scalarFields = {'incmNmse', 'traceRatio', 'largestEigenvalueRatio', ...
    'normalizedEigenvalueProfileError', 'conditionNumber', ...
    'interferenceSubspaceDistance'};
subs = num2cell(indices);
for fieldIndex = 1:numel(scalarFields)
    field = scalarFields{fieldIndex};
    output.(field)(subs{:}) = metrics.(field);
end
if ~isempty(local.peakValue)
    localSubs = [subs, {':'}];
    output.localPeakValue(localSubs{:}) = reshape(local.peakValue, ...
        [ones(1, numel(indices)), numel(local.peakValue)]);
    output.localIntegratedEnergy(localSubs{:}) = reshape(local.integratedEnergy, ...
        [ones(1, numel(indices)), numel(local.integratedEnergy)]);
    output.localPeakWidth(localSubs{:}) = reshape(local.peakWidth, ...
        [ones(1, numel(indices)), numel(local.peakWidth)]);
end
end

function metrics = incmMetrics(estimate, truth, interfererCount)
%INCMMETRICS 比较已形成的INCM；整体尺度与特征结构分别报告。
estimate = (estimate + estimate') / 2;
truth = (truth + truth') / 2;
estimateEigenvalues = sort(max(real(eig(estimate)), 0), 'descend');
truthEigenvalues = sort(max(real(eig(truth)), 0), 'descend');
metrics = struct( ...
    'incmNmse', norm(estimate - truth, 'fro')^2 / norm(truth, 'fro')^2, ...
    'traceRatio', real(trace(estimate)) / real(trace(truth)), ...
    'largestEigenvalueRatio', estimateEigenvalues(1) / truthEigenvalues(1), ...
    'normalizedEigenvalueProfileError', norm( ...
        estimateEigenvalues / sum(estimateEigenvalues) - ...
        truthEigenvalues / sum(truthEigenvalues), 2), ...
    'conditionNumber', cond(estimate), ...
    'interferenceSubspaceDistance', subspaceDistance( ...
        estimate, truth, interfererCount));
end

function distance = subspaceDistance(estimate, truth, interfererCount)
%SUBSPACEDISTANCE J=0时定义为0，J>=M时比较完整空间仍安全返回0。
dimension = min(interfererCount, size(estimate, 1));
if dimension == 0
    distance = 0;
    return;
end
[estimateVectors, estimateValues] = eig(estimate, 'vector');
[~, estimateOrder] = sort(real(estimateValues), 'descend');
[truthVectors, truthValues] = eig(truth, 'vector');
[~, truthOrder] = sort(real(truthValues), 'descend');
estimateBasis = estimateVectors(:, estimateOrder(1:dimension));
truthBasis = truthVectors(:, truthOrder(1:dimension));
distance = norm(estimateBasis * estimateBasis' - ...
    truthBasis * truthBasis', 'fro');
end

function local = localRecDiagnostics(reconstruction, grid, trueNodes, ...
        gridSize, radiusBins)
%LOCALRECDIAGNOSTICS 用真实方向邻域评价原Rec谱，不参与原Rec积分。
included = grid.includedMask;
nodes = grid.nodes(included, :);
powers = reshape(reconstruction.caponPower, [], 1);
energies = reshape(reconstruction.weightedPower, [], 1);
local = localSpectrumDiagnostics(nodes, powers, energies, trueNodes, ...
    gridSize, radiusBins);
end

function local = localCovfitDiagnostics(reconstruction, trueNodes, ...
        gridSize, radiusBins)
%LOCALCOVFITDIAGNOSTICS 用真实方向邻域评价拟合候选，不改变候选集合。
nodes = reconstruction.interfererCandidateNodes;
powers = reshape(reconstruction.interfererFittedPowers, [], 1);
local = localSpectrumDiagnostics(nodes, powers, powers, trueNodes, ...
    gridSize, radiusBins);
end

function local = localSpectrumDiagnostics(nodes, powers, energies, trueNodes, ...
        gridSize, radiusBins)
%LOCALSPECTRUMDIAGNOSTICS 报告邻域峰值、能量及功率加权RMS峰宽。
interfererCount = size(trueNodes, 1);
local = struct('peakValue', nan(1, interfererCount), ...
    'integratedEnergy', nan(1, interfererCount), ...
    'peakWidth', nan(1, interfererCount));
% 邻域半径复用算法配置的NMS半径，并加半个单元覆盖真值落在格点间的情况。
radius = 2 * (radiusBins + 0.5) / gridSize;
for interfererIndex = 1:interfererCount
    distances = vecnorm(nodes - trueNodes(interfererIndex, :), 2, 2);
    mask = distances <= radius & isfinite(powers) & isfinite(energies);
    if ~any(mask)
        continue;
    end
    local.peakValue(interfererIndex) = max(powers(mask));
    local.integratedEnergy(interfererIndex) = sum(energies(mask));
    positivePower = max(real(powers(mask)), 0);
    if sum(positivePower) > 0
        local.peakWidth(interfererIndex) = sqrt(sum(positivePower .* ...
            distances(mask).^2) / sum(positivePower));
    end
end
end

function [desiredOutside, interfererInside] = sectorViolations(obs, cfg)
%SECTORVIOLATIONS 只作实验后评价，不反馈给排除区或候选选择。
nominalDesired = cfg.derived.nominalArrayUnitVectors(:, 1);
actualDirections = obs.truth.arrayUnitVectors;
cosines = max(-1, min(1, real(nominalDesired' * actualDirections)));
distances = acos(cosines);
halfWidth = deg2rad(cfg.spectrum.desiredExclusionHalfWidthDeg);
desiredOutside = distances(1) > halfWidth;
interfererInside = distances(2:end) <= halfWidth;
end

function summary = summarizeDiagnosis(zeroRaw, jitterRaw, cfg, diagnosis, ...
        snapshots, grids, measures, jitterPoints, methods)
%SUMMARIZEDIAGNOSIS 先在线性SINR域统计，再把统计量统一转换为dB。
zeroGroups = numel(snapshots) * numel(grids) * numel(measures) * numel(methods);
zeroSamples = reshape(zeroRaw.sinrLinear, diagnosis.monteCarloTrials, zeroGroups);
zeroSubstream = uav.experiment.substream(cfg, diagnosis.experimentIndex, ...
    1, 0, 0, "bootstrap");
summary.zeroJitter.sinr = reshapeLinearSummary( ...
    uav.experiment.summarizeLinear(zeroSamples, cfg, zeroSubstream), ...
    [numel(snapshots), numel(grids), numel(measures), numel(methods)]);
summary.zeroJitter.originalRec = summarizeIncmDiagnostics( ...
    zeroRaw.originalRec, diagnosis.monteCarloTrials, ...
    [numel(snapshots), numel(grids), numel(measures)], cfg);
summary.zeroJitter.upaCovfit = summarizeIncmDiagnostics( ...
    zeroRaw.upaCovfit, diagnosis.monteCarloTrials, ...
    [numel(snapshots), numel(grids), numel(measures)], cfg);
summary.zeroJitter.desiredOutsideRate = reshape(mean( ...
    zeroRaw.desiredOutsideExclusion, 1), numel(grids), numel(measures));
summary.zeroJitter.maximumInterfererInsideRate = maximumInterfererRate( ...
    zeroRaw.interfererInsideExclusion, 1, ...
    [numel(grids), numel(measures)]);

jitterGroups = numel(jitterPoints) * numel(snapshots) * numel(methods);
jitterSamples = reshape(jitterRaw.sinrLinear, diagnosis.monteCarloTrials, ...
    diagnosis.periods, jitterGroups);
jitterSubstream = uav.experiment.substream(cfg, diagnosis.experimentIndex, ...
    2, 0, 0, "bootstrap");
summary.jitterSweep.sinr = reshapeLinearSummary( ...
    uav.experiment.summarizeTrials(jitterSamples, cfg, jitterSubstream), ...
    [numel(jitterPoints), numel(snapshots), numel(methods)]);
wcoSamples = reshape(jitterRaw.upaCovfitWcoSinrLinear, ...
    diagnosis.monteCarloTrials, diagnosis.periods, ...
    numel(jitterPoints) * numel(snapshots));
summary.jitterSweep.upaCovfitWco = reshapeLinearSummary( ...
    uav.experiment.summarizeTrials(wcoSamples, cfg, jitterSubstream), ...
    [numel(jitterPoints), numel(snapshots)]);
summary.jitterSweep.originalRec = summarizeIncmDiagnostics( ...
    jitterRaw.originalRec, diagnosis.monteCarloTrials * diagnosis.periods, ...
    [numel(jitterPoints), numel(snapshots)], cfg);
summary.jitterSweep.upaCovfit = summarizeIncmDiagnostics( ...
    jitterRaw.upaCovfit, diagnosis.monteCarloTrials * diagnosis.periods, ...
    [numel(jitterPoints), numel(snapshots)], cfg);
summary.jitterSweep.desiredOutsideRate = reshape(mean( ...
    jitterRaw.desiredOutsideExclusion, [1, 2]), numel(jitterPoints), 1);
summary.jitterSweep.maximumInterfererInsideRate = maximumInterfererRate( ...
    jitterRaw.interfererInsideExclusion, [1, 2], [numel(jitterPoints), 1]);
end

function output = reshapeLinearSummary(input, outputShape)
%RESHAPELINEARSUMMARY 保留线性域统计合同并恢复扫描维度。
output = input;
ordinaryFields = {'medianLinear', 'lowLinear', 'meanLinear', ...
    'medianDb', 'lowDb', 'meanDb', 'validCounts'};
for fieldIndex = 1:numel(ordinaryFields)
    field = ordinaryFields{fieldIndex};
    output.(field) = reshape(input.(field), outputShape);
end
output.ci95Linear = reshape(input.ci95Linear, [2, outputShape]);
output.ci95Db = reshape(input.ci95Db, [2, outputShape]);
end

function output = summarizeIncmDiagnostics(raw, rowCount, outputShape, cfg)
%SUMMARIZEINCMDIAGNOSTICS 对NaN不适用值逐组只统计有限样本。
fields = {'incmNmse', 'traceRatio', 'largestEigenvalueRatio', ...
    'normalizedEigenvalueProfileError', 'conditionNumber', ...
    'interferenceSubspaceDistance'};
output = struct();
for fieldIndex = 1:numel(fields)
    field = fields{fieldIndex};
    samples = reshape(raw.(field), rowCount, []);
    finiteSummary = uav.experiment.summarizeFinite(samples, cfg);
    finiteSummary.median = reshape(finiteSummary.median, outputShape);
    finiteSummary.low = reshape(finiteSummary.low, outputShape);
    finiteSummary.mean = reshape(finiteSummary.mean, outputShape);
    finiteSummary.p95 = reshape(finiteSummary.p95, outputShape);
    finiteSummary.validCounts = reshape(finiteSummary.validCounts, outputShape);
    output.(field) = finiteSummary;
end
end

function rate = maximumInterfererRate(values, reductionDimensions, outputShape)
%MAXIMUMINTERFERERRATE 汇总任一真实干扰进入期望排除区的最高发生率。
if size(values, ndims(values)) == 0
    rate = zeros(outputShape);
    return;
end
averaged = mean(values, reductionDimensions);
rate = reshape(max(averaged, [], ndims(averaged)), outputShape);
end

function fig = plotDiagnosis(summary, cfg, snapshots, grids, measures, ...
        jitterPoints, methods, diagnosis)
%PLOTDIAGNOSIS 用五个中文面板展示链路诊断，不据smoke图作性能结论。
fig = uav.plot.newFigure(cfg);
position = fig.Position;
position(3:4) = max(position(3:4), [1600, 900]);
fig.Position = position;
layout = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
labels = methodLabels(methods);
selectedGrid = find(grids == diagnosis.selectedGridSize, 1);
selectedMeasure = find(measures == diagnosis.selectedMeasureName, 1);
colors = lines(numel(methods));

nexttile(layout); hold on;
zeroValues = squeeze(summary.zeroJitter.sinr.meanDb(:, selectedGrid, ...
    selectedMeasure, :));
bar(1:numel(methods), zeroValues.');
set(gca, 'XTick', 1:numel(methods), 'XTickLabel', labels, ...
    'XTickLabelRotation', 20);
xlabel('方法'); ylabel('平均输出SINR/dB');
title('零抖动双快拍对比');
legend("L=" + string(snapshots), 'Location', 'best'); grid on; box on;

nexttile(layout); hold on;
convergenceMethods = [1, 4, 5];
for methodIndex = convergenceMethods
    values = squeeze(summary.zeroJitter.sinr.meanDb(end, :, ...
        selectedMeasure, methodIndex));
    plot(grids, values, '-o', 'LineWidth', cfg.plot.lineWidth, ...
        'MarkerSize', cfg.plot.markerSize, 'Color', colors(methodIndex, :));
end
xlabel('二维谱网格边长'); ylabel('平均输出SINR/dB');
title("L=" + snapshots(end) + "网格收敛");
legend(labels(convergenceMethods), 'Location', 'best'); grid on; box on;

nexttile(layout); hold on;
recError = squeeze(summary.zeroJitter.originalRec. ...
    normalizedEigenvalueProfileError.mean(end, :, selectedMeasure));
upaError = squeeze(summary.zeroJitter.upaCovfit. ...
    normalizedEigenvalueProfileError.mean(end, :, selectedMeasure));
plot(grids, recError, '-o', 'LineWidth', cfg.plot.lineWidth, ...
    'MarkerSize', cfg.plot.markerSize);
plot(grids, upaError, '-s', 'LineWidth', cfg.plot.lineWidth, ...
    'MarkerSize', cfg.plot.markerSize);
xlabel('二维谱网格边长'); ylabel('归一化特征值谱误差');
title('INCM相对谱结构诊断');
legend(["原始Rec", "UPA协方差拟合"], 'Location', 'best'); grid on; box on;

nexttile(layout); hold on;
for methodIndex = 1:numel(methods)
    plot(jitterPoints, squeeze(summary.jitterSweep.sinr.meanDb(:, end, ...
        methodIndex)), '-o', 'LineWidth', cfg.plot.lineWidth, ...
        'MarkerSize', cfg.plot.markerSize, 'Color', colors(methodIndex, :));
end
plot(jitterPoints, squeeze(summary.jitterSweep.upaCovfitWco.meanDb(:, end)), ...
    '--d', 'LineWidth', cfg.plot.lineWidth, 'MarkerSize', cfg.plot.markerSize, ...
    'Color', [0.2, 0.2, 0.2]);
xlabel('三轴姿态抖动RMS/°'); ylabel('平均输出SINR/dB');
title("L=" + snapshots(end) + "抖动扫描");
legend([labels, "UPA协方差拟合+WCO权值"], 'Location', 'best'); grid on; box on;

nexttile(layout); hold on;
plot(jitterPoints, summary.jitterSweep.desiredOutsideRate, '-o', ...
    'LineWidth', cfg.plot.lineWidth, 'MarkerSize', cfg.plot.markerSize);
plot(jitterPoints, summary.jitterSweep.maximumInterfererInsideRate, '-s', ...
    'LineWidth', cfg.plot.lineWidth, 'MarkerSize', cfg.plot.markerSize);
xlabel('三轴姿态抖动RMS/°'); ylabel('扇区违约率');
title('固定期望排除区后验评价');
legend(["期望移出", "任一干扰进入"], 'Location', 'best'); grid on; box on;

allAxes = findall(fig, 'Type', 'axes');
set(allAxes, 'FontSize', cfg.plot.fontSize, 'LineWidth', ...
    max(0.5, cfg.plot.lineWidth - 1));
end

function labels = methodLabels(methods)
%METHODLABELS 返回固定五方法的中文图例。
known = ["oracle", "esp", "wco", "gu_leshem2012_rec", "upa_covfit_incm"];
labelsKnown = ["理想先验MVDR", "特征空间法", "最坏情况法", ...
    "Gu–Leshem原始Rec", "UPA协方差拟合INCM"];
[found, indices] = ismember(methods, known);
if ~all(found)
    error('uav:experiment:UnknownUpaDiagnosisLabel', ...
        '实验9的每个方法都必须具有中文图例。');
end
labels = labelsKnown(indices);
end

function output = diagnosisTable(summary, snapshots, grids, measures, ...
        jitterPoints, methods)
%DIAGNOSISTABLE 生成含中文字段说明的零抖动与抖动扫描长表。
zeroRows = numel(snapshots) * numel(grids) * numel(measures) * numel(methods);
jitterRows = numel(jitterPoints) * numel(snapshots) * (numel(methods) + 1);
rowCount = zeroRows + jitterRows;
stage = strings(rowCount, 1);
snapshotCount = nan(rowCount, 1);
gridSize = nan(rowCount, 1);
measureName = strings(rowCount, 1);
jitterRmsDeg = nan(rowCount, 1);
method = strings(rowCount, 1);
meanSinrDb = nan(rowCount, 1);
finiteSampleCount = zeros(rowCount, 1);
row = 0;
for snapshotIndex = 1:numel(snapshots)
    for gridIndex = 1:numel(grids)
        for measureIndex = 1:numel(measures)
            for methodIndex = 1:numel(methods)
                row = row + 1;
                stage(row) = "零抖动网格测度扫描";
                snapshotCount(row) = snapshots(snapshotIndex);
                gridSize(row) = grids(gridIndex);
                measureName(row) = measures(measureIndex);
                jitterRmsDeg(row) = 0;
                method(row) = methods(methodIndex);
                meanSinrDb(row) = summary.zeroJitter.sinr.meanDb( ...
                    snapshotIndex, gridIndex, measureIndex, methodIndex);
                finiteSampleCount(row) = summary.zeroJitter.sinr.validCounts( ...
                    snapshotIndex, gridIndex, measureIndex, methodIndex);
            end
        end
    end
end
for jitterIndex = 1:numel(jitterPoints)
    for snapshotIndex = 1:numel(snapshots)
        for methodIndex = 1:numel(methods) + 1
            row = row + 1;
            stage(row) = "抖动扫描";
            snapshotCount(row) = snapshots(snapshotIndex);
            jitterRmsDeg(row) = jitterPoints(jitterIndex);
            if methodIndex <= numel(methods)
                method(row) = methods(methodIndex);
                meanSinrDb(row) = summary.jitterSweep.sinr.meanDb( ...
                    jitterIndex, snapshotIndex, methodIndex);
                finiteSampleCount(row) = summary.jitterSweep.sinr.validCounts( ...
                    jitterIndex, snapshotIndex, methodIndex);
            else
                method(row) = "upa_covfit_incm_wco_weight";
                meanSinrDb(row) = summary.jitterSweep.upaCovfitWco.meanDb( ...
                    jitterIndex, snapshotIndex);
                finiteSampleCount(row) = ...
                    summary.jitterSweep.upaCovfitWco.validCounts( ...
                    jitterIndex, snapshotIndex);
            end
        end
    end
end
output = table(stage, snapshotCount, gridSize, measureName, jitterRmsDeg, ...
    method, meanSinrDb, finiteSampleCount);
output.Properties.VariableDescriptions = { ...
    '实验分支', '每周期快拍数', '二维谱网格边长', '积分测度名称', ...
    '三轴姿态抖动RMS（度）', '方法内部名称', ...
    '在线性SINR样本上求均值后转换的dB值', '有限SINR样本数'};
end
