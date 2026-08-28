function [alignedPrevious, diagnostics] = align2D(currentMap, previousMap, mask, cfg)
%ALIGN2D Align a previous spectrum onto current grid coordinates.

validateInputs(currentMap, previousMap, mask, cfg);
maxShift = reshape(cfg.spectrum.alignmentMaxShiftBins, 1, 2);
shifts = zeros(prod(2 * maxShift + 1), 2);
scores = -inf(size(shifts, 1), 1);
shiftIndex = 0;
for rowShift = -maxShift(1):maxShift(1)
    for columnShift = -maxShift(2):maxShift(2)
        shiftIndex = shiftIndex + 1;
        shifts(shiftIndex, :) = [rowShift, columnShift];
        [shifted, shiftedMask] = shiftWithoutWrap(previousMap, mask, ...
            rowShift, columnShift);
        overlap = mask & shiftedMask;
        scores(shiftIndex) = normalizedCorrelation(currentMap(overlap), ...
            shifted(overlap), cfg.spectrum.alignmentMinimumOverlap);
    end
end

[bestScore, order] = sort(scores, 'descend');
bestIndex = order(1);
if numel(bestScore) > 1
    secondScore = bestScore(2);
else
    secondScore = -inf;
end
weakEvidence = ~isfinite(bestScore(1)) || ...
    bestScore(1) < cfg.spectrum.alignmentMinCorrelation;
ambiguousEvidence = isfinite(secondScore) && ...
    bestScore(1) - secondScore < cfg.spectrum.alignmentAmbiguityMargin;
if weakEvidence || ambiguousEvidence
    selectedShift = [0, 0];
    alignedPrevious = previousMap;
    if weakEvidence
        fallbackReason = "weak-evidence";
    else
        fallbackReason = "ambiguous-evidence";
    end
    usedIdentity = true;
else
    selectedShift = shifts(bestIndex, :);
    alignedPrevious = shiftWithoutWrap(previousMap, true(size(mask)), ...
        selectedShift(1), selectedShift(2));
    fallbackReason = "none";
    usedIdentity = all(selectedShift == 0);
end
alignedPrevious(~mask) = 0;
diagnostics = struct( ...
    'shiftBins', selectedShift, ...
    'bestCorrelation', bestScore(1), ...
    'secondCorrelation', secondScore, ...
    'usedIdentity', usedIdentity, ...
    'fallbackReason', fallbackReason, ...
    'candidateShifts', shifts, ...
    'candidateCorrelations', scores);
end

function validateInputs(currentMap, previousMap, mask, cfg)
if ~ismatrix(currentMap) || ~isequal(size(currentMap), size(previousMap), size(mask)) || ...
        ~islogical(mask) || any(~isfinite(currentMap(mask))) || ...
        any(~isfinite(previousMap(mask))) || any(currentMap(mask) < 0) || ...
        any(previousMap(mask) < 0)
    error('uav:spectrum:InvalidAlignmentInput', ...
        'Spectra and logical mask must have equal sizes and finite nonnegative support.');
end
if ~isstruct(cfg) || ~isfield(cfg, 'spectrum') || ...
        ~isfield(cfg.spectrum, 'alignmentMaxShiftBins')
    error('uav:spectrum:InvalidConfiguration', ...
        'cfg is missing spectrum alignment fields.');
end
end

function score = normalizedCorrelation(first, second, minimumOverlap)
if numel(first) < minimumOverlap
    score = -inf;
    return;
end
first = real(first(:));
second = real(second(:));
first = first - mean(first);
second = second - mean(second);
denominator = norm(first) * norm(second);
if denominator <= eps(max(1, norm(first) + norm(second)))
    score = -inf;
else
    score = max(-1, min(1, (first' * second) / denominator));
end
end

function [shifted, shiftedMask] = shiftWithoutWrap(values, mask, rowShift, columnShift)
[rowCount, columnCount] = size(values);
shifted = zeros(size(values), 'like', values);
shiftedMask = false(size(mask));
sourceRows = max(1, 1 - rowShift):min(rowCount, rowCount - rowShift);
sourceColumns = max(1, 1 - columnShift):min(columnCount, columnCount - columnShift);
if isempty(sourceRows) || isempty(sourceColumns)
    return;
end
targetRows = sourceRows + rowShift;
targetColumns = sourceColumns + columnShift;
shifted(targetRows, targetColumns) = values(sourceRows, sourceColumns);
shiftedMask(targetRows, targetColumns) = mask(sourceRows, sourceColumns);
end
