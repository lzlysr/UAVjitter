function summary = summarizeFinite(values, cfg)
%SUMMARIZEFINITE Aggregate finite values while preserving unavailable NaNs.

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
