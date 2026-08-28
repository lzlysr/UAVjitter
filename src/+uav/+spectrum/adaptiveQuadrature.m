function quadrature = adaptiveQuadrature(powerMap, spectrumDiagnostics, includedMask, peakIndices, cfg)
%ADAPTIVEQUADRATURE 对峰值邻域细分、其余干扰区域粗聚合。
%   输入功率图和掩码与Capon网格同尺寸；peakIndices是selectPeaks给出的
%   二维线性索引。输出quadrature包含N×2节点、N×1面积权重和谱值。
%
%   设计目标是在峰值附近保留角域分辨率，同时减少平坦背景区域的节点数。
%   细区和粗区由互斥掩码划分，因此不会有意重复表示同一原始网格单元。

validateInputs(powerMap, spectrumDiagnostics, includedMask, peakIndices, cfg);
gridStep = sqrt(spectrumDiagnostics.cellArea);
fineMask = false(size(includedMask));
radius = cfg.spectrum.quadratureFineRadiusBins;
% 把每个检测峰周围指定半径的网格单元标为细分区域。
for peakIndex = reshape(peakIndices, 1, [])
    [peakRow, peakColumn] = ind2sub(size(includedMask), peakIndex);
    rows = max(1, peakRow - radius):min(size(includedMask, 1), peakRow + radius);
    columns = max(1, peakColumn - radius):min(size(includedMask, 2), peakColumn + radius);
    fineMask(rows, columns) = true;
end
fineMask = fineMask & includedMask;

% 峰值单元内部生成refinement×refinement子节点；非峰值区按stride聚块。
[fineNodes, fineWeights, fineValues] = refineCells(powerMap, ...
    spectrumDiagnostics, fineMask, gridStep, ...
    cfg.spectrum.quadratureRefinement, cfg);
[coarseNodes, coarseWeights, coarseValues] = aggregateCoarseCells(powerMap, ...
    spectrumDiagnostics, includedMask & ~fineMask, spectrumDiagnostics.cellArea, ...
    cfg.spectrum.quadratureCoarseStride);

nodes = [fineNodes; coarseNodes];
weights = [fineWeights; coarseWeights];
values = [fineValues; coarseValues];
[nodes, weights, values, mergeDiagnostics] = ...
    uav.spectrum.mergeQuadratureNodes( ...
    nodes, weights, values, ...
    cfg.spectrum.quadratureDedupTolerance);
quadrature = struct( ...
    'nodes', nodes, ...
    'weights', weights, ...
    'values', values, ...
    'diagnostics', struct( ...
        'mode', "adaptive", ...
        'fineCellCount', nnz(fineMask), ...
        'coarseCellCount', nnz(includedMask & ~fineMask), ...
        'fineNodeCount', size(fineNodes, 1), ...
        'coarseNodeCount', size(coarseNodes, 1), ...
        'deduplicatedNodeCount', mergeDiagnostics.deduplicatedNodeCount, ...
        'representedArea', sum(weights), ...
        'fineMask', fineMask));
end

function [nodes, weights, values] = refineCells( ...
        powerMap, diagnostics, mask, gridStep, refinement, cfg)
linearIndices = find(mask);
nodes = zeros(0, 2);
weights = zeros(0, 1);
values = zeros(0, 1);
offsets = (((1:refinement) - 0.5) / refinement - 0.5) * gridStep;
for linearIndex = reshape(linearIndices, 1, [])
    % 子节点位于原单元内部的等距中心位置，而不是落在单元边界上。
    center = [diagnostics.xiX(linearIndex), diagnostics.xiY(linearIndex)];
    [offsetX, offsetY] = meshgrid(offsets, offsets);
    cellNodes = [center(1) + offsetX(:), center(2) + offsetY(:)];
    admissible = uav.spectrum.admissibleDirections(cellNodes, cfg);
    cellNodes = cellNodes(admissible, :);
    if isempty(cellNodes)
        % 极端边界单元若所有子点不可接受，保留已验证可接受的原中心。
        cellNodes = center;
    end
    cellWeights = diagnostics.cellArea / size(cellNodes, 1) * ...
        ones(size(cellNodes, 1), 1);
    cellValues = interp2(diagnostics.xiX, diagnostics.xiY, powerMap, ...
        cellNodes(:, 1), cellNodes(:, 2), 'linear');
    cellValues(~isfinite(cellValues)) = powerMap(linearIndex);
    nodes = [nodes; cellNodes]; %#ok<AGROW>
    weights = [weights; cellWeights]; %#ok<AGROW>
    values = [values; max(real(cellValues), 0)]; %#ok<AGROW>
end
end

function [nodes, weights, values] = aggregateCoarseCells(powerMap, diagnostics, mask, cellArea, stride)
[rows, columns] = find(mask);
if isempty(rows)
    nodes = zeros(0, 2);
    weights = zeros(0, 1);
    values = zeros(0, 1);
    return;
end
blockRows = floor((rows - 1) / stride);
blockColumns = floor((columns - 1) / stride);
% 相同粗块内的原始单元合并成一个代表节点，总面积保持不变。
[~, ~, groups] = unique([blockRows, blockColumns], 'rows', 'stable');
groupCount = max(groups);
nodes = zeros(groupCount, 2);
weights = zeros(groupCount, 1);
values = zeros(groupCount, 1);
linearIndices = sub2ind(size(mask), rows, columns);
for group = 1:groupCount
    members = groups == group;
    count = nnz(members);
    memberIndices = linearIndices(members);
    memberNodes = [diagnostics.xiX(memberIndices), diagnostics.xiY(memberIndices)];
    centroid = mean(memberNodes, 1);
    % 选择离几何中心最近的真实网格点，避免代表点落到不可见域外。
    [~, representativeIndex] = min(sum((memberNodes - centroid).^2, 2));
    nodes(group, :) = memberNodes(representativeIndex, :);
    weights(group) = count * cellArea;
    values(group) = mean(real(powerMap(linearIndices(members))));
end
end

function validateInputs(powerMap, diagnostics, includedMask, peakIndices, cfg)
required = {'xiX', 'xiY', 'visibleMask', 'cellArea'};
if ~ismatrix(powerMap) || ~isstruct(diagnostics) || ~all(isfield(diagnostics, required)) || ...
        ~isequal(size(powerMap), size(includedMask), size(diagnostics.xiX), ...
            size(diagnostics.xiY), size(diagnostics.visibleMask)) || ...
        ~islogical(includedMask) || any(~isfinite(powerMap(includedMask))) || ...
        any(powerMap(includedMask) < 0) || any(includedMask & ~diagnostics.visibleMask, 'all') || ...
        ~isvector(peakIndices) || any(peakIndices < 1) || ...
        any(peakIndices ~= floor(peakIndices)) || any(peakIndices > numel(powerMap)) || ...
        any(~includedMask(peakIndices))
    error('uav:spectrum:InvalidQuadratureInput', ...
        'Quadrature inputs must describe finite included visible grid cells.');
end
if ~isstruct(cfg) || ~isfield(cfg, 'spectrum')
    error('uav:spectrum:InvalidConfiguration', 'cfg is missing quadrature fields.');
end
gridNodes = [diagnostics.xiX(includedMask), diagnostics.xiY(includedMask)];
if ~all(uav.spectrum.admissibleDirections(gridNodes, cfg))
    error('uav:spectrum:InvalidQuadratureInput', ...
        'Every included cell center must satisfy continuous admissibility.');
end
end
