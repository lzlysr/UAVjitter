function summary = summarizeLinear(values, cfg, bootstrapSubstream)
%SUMMARIZELINEAR 在线性功率域汇总正值样本并进行配对bootstrap。
%   values为[样本行×方法]。每次重采样对所有方法使用同一组行索引，保持
%   方法间配对关系。该函数适用于独立样本行；含多周期相关性的主实验应
%   调用summarizeTrials，以“试次”作为重采样簇。

if ~isnumeric(values) || ~ismatrix(values) || isempty(values) || ...
        any(values(isfinite(values)) <= 0)
    error('uav:experiment:InvalidLinearSamples', ...
        'Linear samples must be a nonempty matrix of positive values or NaN.');
end
if any(~isfinite(values) & ~isnan(values), 'all')
    error('uav:experiment:InvalidLinearSamples', ...
        'Linear samples may use NaN only for unavailable values.');
end
groupCount = size(values, 2);
medianLinear = nan(1, groupCount);
lowLinear = nan(1, groupCount);
meanLinear = nan(1, groupCount);
ci95Linear = nan(2, groupCount);
validCounts = zeros(1, groupCount);
stream = uav.util.makeStream(cfg.random.baseSeed, bootstrapSubstream);
rowCount = size(values, 1);
bootstrapMedians = nan(cfg.experiments.bootstrapReplicates, groupCount);
for replicate = 1:cfg.experiments.bootstrapReplicates
    % 每次只生成一列索引，避免paper规模构造巨大的三维索引数组。
    sharedIndices = randi(stream, rowCount, rowCount, 1);
    for group = 1:groupCount
        selectedValues = values(sharedIndices, group);
        selectedValues = selectedValues(isfinite(selectedValues));
        if ~isempty(selectedValues)
            bootstrapMedians(replicate, group) = median(selectedValues);
        end
    end
end
for group = 1:groupCount
    finiteValues = values(isfinite(values(:, group)), group);
    validCounts(group) = numel(finiteValues);
    if isempty(finiteValues)
        continue;
    end
    medianLinear(group) = quantile(finiteValues, ...
        cfg.experiments.aggregation.medianQuantile);
    lowLinear(group) = quantile(finiteValues, ...
        cfg.experiments.aggregation.lowQuantile);
    meanLinear(group) = sum(finiteValues) / numel(finiteValues);
    finiteBootstrap = bootstrapMedians(isfinite(bootstrapMedians(:, group)), group);
    ci95Linear(:, group) = quantile(finiteBootstrap, ...
        cfg.experiments.aggregation.confidenceQuantiles).';
end
summary = struct( ...
    'medianLinear', medianLinear, ...
    'lowLinear', lowLinear, ...
    'meanLinear', meanLinear, ...
    'ci95Linear', ci95Linear, ...
    'medianDb', linearToDb(medianLinear), ...
    'lowDb', linearToDb(lowLinear), ...
    'meanDb', linearToDb(meanLinear), ...
    'ci95Db', linearToDb(ci95Linear), ...
    'validCounts', validCounts, ...
    'bootstrapSubstream', bootstrapSubstream, ...
    'bootstrapReplicates', cfg.experiments.bootstrapReplicates, ...
    'bootstrapUnit', "row", ...
    'aggregationDomain', "linear");
end

function valuesDb = linearToDb(values)
% 先在线性域完成中位数/均值/置信区间，再把最终统计量转换为dB。
valuesDb = nan(size(values));
finite = isfinite(values);
valuesDb(finite) = 10 * log10(values(finite));
end
