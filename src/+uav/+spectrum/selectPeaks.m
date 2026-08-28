function [peaks, diagnostics] = selectPeaks(powerMap, spectrumDiagnostics, includedMask, cfg)
%SELECTPEAKS Select spectrum peaks by thresholded nonmaximum suppression.

validateInputs(powerMap, spectrumDiagnostics, includedMask, cfg);
cap = max(cfg.derived.M - 1, 0);
eligibleIndices = find(includedMask);
eligiblePowers = real(powerMap(eligibleIndices));
if isempty(eligiblePowers) || cap == 0 || max(eligiblePowers) <= 0
    selectedIndices = zeros(0, 1);
    threshold = 0;
else
    threshold = cfg.spectrum.peakThresholdRatio * max(eligiblePowers);
    keep = eligiblePowers >= threshold;
    candidates = eligibleIndices(keep);
    candidatePowers = eligiblePowers(keep);
    [~, order] = sort(candidatePowers, 'descend');
    candidates = candidates(order);
    selectedIndices = zeros(0, 1);
    radius = cfg.spectrum.peakNmsRadiusBins;
    [candidateRows, candidateColumns] = ind2sub(size(powerMap), candidates);
    for index = 1:numel(candidates)
        if isempty(selectedIndices)
            isSuppressed = false;
        else
            [selectedRows, selectedColumns] = ind2sub(size(powerMap), selectedIndices);
            isSuppressed = any(max(abs(selectedRows - candidateRows(index)), ...
                abs(selectedColumns - candidateColumns(index))) <= radius);
        end
        if ~isSuppressed
            selectedIndices(end + 1, 1) = candidates(index); %#ok<AGROW>
            if numel(selectedIndices) == cap
                break;
            end
        end
    end
end

peaks = [spectrumDiagnostics.xiX(selectedIndices), ...
    spectrumDiagnostics.xiY(selectedIndices)];
diagnostics = struct( ...
    'linearIndices', selectedIndices, ...
    'powers', powerMap(selectedIndices), ...
    'threshold', threshold, ...
    'cap', cap, ...
    'candidateCount', numel(eligibleIndices));
end

function validateInputs(powerMap, spectrumDiagnostics, includedMask, cfg)
required = {'xiX', 'xiY', 'visibleMask'};
if ~ismatrix(powerMap) || ~isstruct(spectrumDiagnostics) || ...
        ~all(isfield(spectrumDiagnostics, required)) || ...
        ~isequal(size(powerMap), size(includedMask), size(spectrumDiagnostics.xiX), ...
            size(spectrumDiagnostics.xiY), size(spectrumDiagnostics.visibleMask)) || ...
        ~islogical(includedMask) || any(~isfinite(powerMap(includedMask))) || ...
        any(powerMap(includedMask) < 0) || any(includedMask & ~spectrumDiagnostics.visibleMask, 'all')
    error('uav:spectrum:InvalidPeakInput', ...
        'Peak input maps and masks must be finite, nonnegative, visible, and size matched.');
end
if ~isstruct(cfg) || ~isfield(cfg, 'derived') || ~isfield(cfg, 'spectrum')
    error('uav:spectrum:InvalidConfiguration', 'cfg is missing peak-selection fields.');
end
end
