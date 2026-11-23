function summary = summarizePoints(values, cfg, experimentIndex, comparedMethods)
%SUMMARIZEPOINTS 对每个扫描点独立执行试次簇统计。
%   输入维度[试次×周期×扫描点×方法]；每个扫描点使用独立bootstrap
%   子流，保证结果可复现且一个点的重采样不会影响另一个点。

dimensions = size(values);
if numel(dimensions) < 4
    dimensions(end + 1:4) = 1;
end
pointCount = dimensions(3);
groupCount = dimensions(4);
summary = struct( ...
    'medianLinear', nan(pointCount, groupCount), ...
    'lowLinear', nan(pointCount, groupCount), ...
    'meanLinear', nan(pointCount, groupCount), ...
    'ci95Linear', nan(2, pointCount, groupCount), ...
    'medianDb', nan(pointCount, groupCount), ...
    'lowDb', nan(pointCount, groupCount), ...
    'meanDb', nan(pointCount, groupCount), ...
    'ci95Db', nan(2, pointCount, groupCount), ...
    'validCounts', zeros(pointCount, groupCount), ...
    'aggregationDomain', "linear");
if nargin < 4, comparedMethods = []; end
if ~isempty(comparedMethods)
    summary.pairedDifferenceDb = nan(pointCount,1);
    summary.pairedDifferenceCi95Db = nan(pointCount,2);
end
for pointIndex = 1:pointCount
    samples = reshape(values(:, :, pointIndex, :), ...
        dimensions(1), dimensions(2), groupCount);
    bootstrapSubstream = uav.experiment.substream(cfg, experimentIndex, ...
        pointIndex, 0, 0, "bootstrap");
    pointSummary = uav.experiment.summarizeTrials( ...
        samples, cfg, bootstrapSubstream, comparedMethods);
    summary.medianLinear(pointIndex, :) = pointSummary.medianLinear;
    summary.lowLinear(pointIndex, :) = pointSummary.lowLinear;
    summary.meanLinear(pointIndex, :) = pointSummary.meanLinear;
    summary.ci95Linear(:, pointIndex, :) = ...
        reshape(pointSummary.ci95Linear, 2, 1, groupCount);
    summary.medianDb(pointIndex, :) = pointSummary.medianDb;
    summary.lowDb(pointIndex, :) = pointSummary.lowDb;
    summary.meanDb(pointIndex, :) = pointSummary.meanDb;
    summary.ci95Db(:, pointIndex, :) = ...
        reshape(pointSummary.ci95Db, 2, 1, groupCount);
    summary.validCounts(pointIndex, :) = pointSummary.validCounts;
    if ~isempty(comparedMethods)
        summary.pairedDifferenceDb(pointIndex) = pointSummary.pairedDifferenceDb;
        summary.pairedDifferenceCi95Db(pointIndex,:) = pointSummary.pairedDifferenceCi95Db.';
    end
end
end
