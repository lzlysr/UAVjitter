function [peaks, diagnostics] = selectPeaks(powerMap, spectrumDiagnostics, includedMask, cfg)
%SELECTPEAKS 对增强谱执行阈值筛选和二维非极大值抑制。
%   峰值数量不读取真实干扰数J，只由相对阈值、抑制半径和M-1上限决定。
%   peaks输出Npeak×2方向余弦，diagnostics保留二维线性索引和峰值功率。

validateInputs(powerMap, spectrumDiagnostics, includedMask, cfg);
cap = max(cfg.derived.M - 1, 0);
eligibleIndices = find(includedMask);
eligiblePowers = real(powerMap(eligibleIndices));
if isempty(eligiblePowers) || cap == 0 || max(eligiblePowers) <= 0
    selectedIndices = zeros(0, 1);
    threshold = 0;
else
    % 只保留达到全局候选最大值一定比例的谱单元。
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
        % Chebyshev网格距离不超过radius的较弱候选被已有强峰抑制。
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
