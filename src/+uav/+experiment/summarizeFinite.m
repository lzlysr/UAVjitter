function summary = summarizeFinite(values, cfg)
%SUMMARIZEFINITE 汇总允许含NaN的普通有限指标。
%   values每列代表一个方法或变体。NaN表示该指标对该方法不可用；Inf和
%   -Inf始终视为错误。输出中位数、低分位、均值、95%运行时间和有效数。

if ~isnumeric(values) || ~ismatrix(values) || isempty(values) || ...
        any(~isfinite(values) & ~isnan(values), 'all')
    error('uav:experiment:InvalidFiniteSamples', ...
        'Samples must be a nonempty numeric matrix with NaN only for unavailable values.');
end
groupCount = size(values, 2);
medianValue = nan(1, groupCount);
lowValue = nan(1, groupCount);
meanValue = nan(1, groupCount);
runtime95 = nan(1, groupCount);
validCounts = zeros(1, groupCount);
for group = 1:groupCount
    % 每组独立剔除NaN，不能用零填充扭曲统计量。
    finiteValues = values(isfinite(values(:, group)), group);
    validCounts(group) = numel(finiteValues);
    if isempty(finiteValues)
        continue;
    end
    medianValue(group) = quantile(finiteValues, ...
        cfg.experiments.aggregation.medianQuantile);
    lowValue(group) = quantile(finiteValues, ...
        cfg.experiments.aggregation.lowQuantile);
    meanValue(group) = sum(finiteValues) / numel(finiteValues);
    runtime95(group) = quantile(finiteValues, ...
        cfg.experiments.aggregation.runtimeQuantile);
end
summary = struct('median', medianValue, 'low', lowValue, ...
    'mean', meanValue, 'p95', runtime95, 'validCounts', validCounts);
end
