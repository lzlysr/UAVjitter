function [mergedNodes, mergedWeights, mergedValues, diagnostics] = ...
        mergeQuadratureNodes(nodes, weights, values, tolerance)
%MERGEQUADRATURENODES Merge near-duplicate nodes without quantizing coordinates.

validateInputs(nodes, weights, values, tolerance);
if isempty(nodes)
    mergedNodes = nodes;
    mergedWeights = weights;
    mergedValues = values;
    diagnostics = struct('deduplicatedNodeCount', 0, ...
        'inputNodeCount', 0, 'outputNodeCount', 0);
    return;
end

keys = round(nodes / tolerance);
[~, firstIndices, groups] = unique(keys, 'rows', 'stable');
mergedWeights = accumarray(groups, weights);
mergedValues = accumarray(groups, weights .* values) ./ mergedWeights;
mergedNodes = nodes(firstIndices, :);
diagnostics = struct( ...
    'deduplicatedNodeCount', size(nodes, 1) - size(mergedNodes, 1), ...
    'inputNodeCount', size(nodes, 1), ...
    'outputNodeCount', size(mergedNodes, 1));
end

function validateInputs(nodes, weights, values, tolerance)
nodeCount = size(nodes, 1);
if ~isnumeric(nodes) || size(nodes, 2) ~= 2 || any(~isfinite(nodes), 'all') || ...
        ~isnumeric(weights) || ~isvector(weights) || numel(weights) ~= nodeCount || ...
        any(~isfinite(weights)) || any(weights <= 0) || ...
        ~isnumeric(values) || ~isvector(values) || numel(values) ~= nodeCount || ...
        any(~isfinite(values)) || any(values < 0) || ...
        ~isscalar(tolerance) || ~isfinite(tolerance) || tolerance <= 0
    error('uav:spectrum:InvalidQuadratureMergeInput', ...
        'Nodes, positive weights, nonnegative values, and tolerance must be valid.');
end
weights = weights(:); %#ok<NASGU>
values = values(:); %#ok<NASGU>
end
