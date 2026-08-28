function [mergedNodes, mergedWeights, mergedValues, diagnostics] = ...
        mergeQuadratureNodes(nodes, weights, values, tolerance)
%MERGEQUADRATURENODES 合并粗细求积中重复或近重复的方向节点。
%   nodes/weights/values分别为N×2节点、单元面积权重和谱值。
%   tolerance只用于构造分组键；输出节点仍保留首次出现的原坐标，不把
%   实际方向量化到容差网格。合并后的谱值采用面积加权平均。

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
% stable保证输出节点顺序与首次出现顺序一致，便于确定性复现。
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
