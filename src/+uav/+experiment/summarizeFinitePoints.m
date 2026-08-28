function summary = summarizeFinitePoints(values, cfg)
%SUMMARIZEFINITEPOINTS Summarize finite trial-period data by point and group.

dimensions = size(values);
if numel(dimensions) < 4
    dimensions(end + 1:4) = 1;
end
pointCount = dimensions(3);
groupCount = dimensions(4);
summary = struct('median', nan(pointCount, groupCount), ...
    'low', nan(pointCount, groupCount), 'mean', nan(pointCount, groupCount), ...
    'p95', nan(pointCount, groupCount), ...
    'validCounts', zeros(pointCount, groupCount));
for pointIndex = 1:pointCount
    samples = reshape(values(:, :, pointIndex, :), [], groupCount);
    pointSummary = uav.experiment.summarizeFinite(samples, cfg);
    summary.median(pointIndex, :) = pointSummary.median;
    summary.low(pointIndex, :) = pointSummary.low;
    summary.mean(pointIndex, :) = pointSummary.mean;
    summary.p95(pointIndex, :) = pointSummary.p95;
    summary.validCounts(pointIndex, :) = pointSummary.validCounts;
end
end
