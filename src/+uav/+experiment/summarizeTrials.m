function summary = summarizeTrials(values, cfg, bootstrapSubstream)
%SUMMARIZETRIALS Summarize trial-period data with paired trial-cluster bootstrap.

if ~isnumeric(values) || isempty(values) || ndims(values) > 3 || ...
        any(values(isfinite(values)) <= 0) || ...
        any(~isfinite(values) & ~isnan(values), 'all')
    error('uav:experiment:InvalidTrialSamples', ...
        'Trial samples must be positive finite values or NaN.');
end
dimensions = size(values);
if numel(dimensions) < 3
    dimensions(3) = 1;
end
trialCount = dimensions(1);
periodCount = dimensions(2);
groupCount = dimensions(3);
values = reshape(values, trialCount, periodCount, groupCount);
summary = uav.experiment.summarizeLinear(reshape(values, [], groupCount), ...
    cfg, bootstrapSubstream);

stream = uav.util.makeStream(cfg.random.baseSeed, bootstrapSubstream);
replicateCount = cfg.experiments.bootstrapReplicates;
bootstrapMedians = nan(replicateCount, groupCount);
for replicate = 1:replicateCount
    selectedTrials = randi(stream, trialCount, trialCount, 1);
    selectedValues = values(selectedTrials, :, :);
    for group = 1:groupCount
        groupValues = reshape(selectedValues(:, :, group), [], 1);
        groupValues = groupValues(isfinite(groupValues));
        if ~isempty(groupValues)
            bootstrapMedians(replicate, group) = median(groupValues);
        end
    end
end
for group = 1:groupCount
    finiteReplicates = bootstrapMedians(isfinite(bootstrapMedians(:, group)), group);
    if ~isempty(finiteReplicates)
        summary.ci95Linear(:, group) = quantile(finiteReplicates, ...
            cfg.experiments.aggregation.confidenceQuantiles).';
    end
end
summary.ci95Db = linearToDb(summary.ci95Linear);
summary.bootstrapUnit = "trial";
summary.independentTrialCount = trialCount;
summary.periodsPerTrial = periodCount;
end

function valuesDb = linearToDb(values)
valuesDb = nan(size(values));
finite = isfinite(values);
valuesDb(finite) = 10 * log10(values(finite));
end
