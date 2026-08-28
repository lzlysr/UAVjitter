function quadrature = adaptiveQuadrature(powerMap, spectrumDiagnostics, includedMask, peakIndices, cfg)
%ADAPTIVEQUADRATURE Partition grid cells into fine peak and coarse regions.

validateInputs(powerMap, spectrumDiagnostics, includedMask, peakIndices, cfg);
gridStep = sqrt(spectrumDiagnostics.cellArea);
fineMask = false(size(includedMask));
radius = cfg.spectrum.quadratureFineRadiusBins;
for peakIndex = reshape(peakIndices, 1, [])
    [peakRow, peakColumn] = ind2sub(size(includedMask), peakIndex);
    rows = max(1, peakRow - radius):min(size(includedMask, 1), peakRow + radius);
    columns = max(1, peakColumn - radius):min(size(includedMask, 2), peakColumn + radius);
    fineMask(rows, columns) = true;
end
fineMask = fineMask & includedMask;

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
    center = [diagnostics.xiX(linearIndex), diagnostics.xiY(linearIndex)];
    [offsetX, offsetY] = meshgrid(offsets, offsets);
    cellNodes = [center(1) + offsetX(:), center(2) + offsetY(:)];
    admissible = uav.spectrum.admissibleDirections(cellNodes, cfg);
    cellNodes = cellNodes(admissible, :);
    if isempty(cellNodes)
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
