function summary = summarizeLinear(values, cfg, bootstrapSubstream)
%SUMMARIZELINEAR Aggregate positive linear-domain samples with paired bootstrap.

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
valuesDb = nan(size(values));
finite = isfinite(values);
valuesDb(finite) = 10 * log10(values(finite));
end
