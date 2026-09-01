function [RinHat, diagnostics] = upaCovfit(obs, cfg)
%UPACOVFIT 以Capon候选和非负协方差拟合重构UPA干扰加噪声协方差。
%   Capon谱只负责定位少量候选方向；候选功率由样本协方差的非负拟合
%   得到。期望区候选参与联合拟合但不进入最终INCM，因此可降低期望
%   信号自消风险。算法只读取白名单观测与配置先验，不读取任何真值。

validateInputs(obs, cfg);
M = cfg.derived.M;
methodCfg = cfg.incm.staticCovfit;
snapshotCount = size(obs.X, 2);

% 用局部配置覆盖CovFit专用谱参数，不改变调用方的全局配置。
spectrumCfg = cfg;
spectrumCfg.spectrum.diagonalLoading = methodCfg.loadingRatio;
spectrumCfg.spectrum.peakThresholdRatio = methodCfg.peakThresholdRatio;
spectrumCfg.spectrum.peakNmsRadiusBins = methodCfg.peakNmsRadiusBins;
[powerMap, spectrumDiagnostics] = uav.spectrum.loadedCapon(obs.Rx, spectrumCfg);

allNodes = [spectrumDiagnostics.xiX(:), spectrumDiagnostics.xiY(:)];
outsideDesired = uav.spectrum.admissibleDirections(allNodes, cfg);
interfererMask = spectrumDiagnostics.visibleMask & ...
    reshape(outsideDesired, size(spectrumDiagnostics.visibleMask));
desiredMask = spectrumDiagnostics.visibleMask & ~interfererMask;
nominalNode = cfg.derived.nominalDirectionCosines(:, 1).';
tolerance = cfg.spectrum.quadratureDedupTolerance;

% 期望字典始终包含名义方向；再按需加入期望区最强有限谱点。
[desiredPeakNode, desiredPeakStatus] = strongestFiniteNode( ...
    powerMap, spectrumDiagnostics, desiredMask);
[desiredNodes, desiredPeakAdded] = appendUniqueNode( ...
    nominalNode, desiredPeakNode, tolerance);
if desiredPeakStatus == "unavailable"
    desiredFallback = "nominal-only";
elseif desiredPeakAdded
    desiredFallback = "strongest-finite";
else
    desiredFallback = "deduplicated-with-nominal";
end

% 干扰区先做阈值筛选和二维非极大值抑制，再施加空间维数上限。
[peakNodes, peakDiagnostics] = uav.spectrum.selectPeaks( ...
    powerMap, spectrumDiagnostics, interfererMask, spectrumCfg);
candidateCap = min(methodCfg.maximumCandidateCount, M - 1);
peakNodes = peakNodes(1:min(size(peakNodes, 1), candidateCap), :);
interfererFallback = "nms";
if candidateCap == 0
    interfererNodes = zeros(0, 2);
    interfererFallback = "cap-zero";
elseif isempty(peakNodes)
    [fallbackNode, fallbackStatus] = strongestFiniteNode( ...
        powerMap, spectrumDiagnostics, interfererMask);
    if fallbackStatus == "unavailable"
        interfererNodes = zeros(0, 2);
        interfererFallback = "unavailable";
    else
        interfererNodes = fallbackNode;
        interfererFallback = "strongest-finite-fallback";
    end
else
    interfererNodes = peakNodes;
end
interfererNodes = uniqueNodes(interfererNodes, desiredNodes, tolerance);

[noisePower, noiseSource, noiseDiagnostics] = estimateNoisePower( ...
    obs.Rx, snapshotCount, cfg);
candidateNodes = [desiredNodes; interfererNodes];
candidateManifold = uav.array.steeringUPA(candidateNodes.', ...
    cfg.array.Nx, cfg.array.Ny, cfg.array.spacingOverLambda);
[fittedPowers, fitDiagnostics] = uav.incm.fitNonnegativePowers( ...
    obs.Rx, candidateManifold, noisePower, methodCfg.fitRegularization);

% 期望候选只用于分离期望能量；最终INCM仅保留干扰候选和白噪声。
desiredCount = size(desiredNodes, 1);
interfererPowers = fittedPowers(desiredCount + 1:end);
interfererManifold = candidateManifold(:, desiredCount + 1:end);
if isempty(interfererPowers)
    interferenceCovariance = zeros(M);
else
    interferenceCovariance = (interfererManifold .* interfererPowers.') * ...
        interfererManifold';
end
RinHat = interferenceCovariance + noisePower * eye(M);
RinHat = (RinHat + RinHat') / 2;

diagnostics = struct( ...
    'estimator', "upa-covfit", ...
    'loading', spectrumDiagnostics.loading, ...
    'noiseEstimate', noisePower, ...
    'noiseSource', noiseSource, ...
    'noisePower', noisePower, ...
    'effectiveSnapshotCount', snapshotCount, ...
    'noiseDiagnostics', noiseDiagnostics, ...
    'desiredCandidateNodes', desiredNodes, ...
    'interfererCandidateNodes', interfererNodes, ...
    'interfererCandidateCap', candidateCap, ...
    'interfererCandidateCount', size(interfererNodes, 1), ...
    'candidateManifold', candidateManifold, ...
    'fittedPowers', fittedPowers, ...
    'desiredFittedPowers', fittedPowers(1:desiredCount), ...
    'interfererFittedPowers', interfererPowers, ...
    'fitResidual', fitDiagnostics.relativeResidual, ...
    'fitDiagnostics', fitDiagnostics, ...
    'conditionNumber', cond(RinHat), ...
    'desiredCandidateFallback', desiredFallback, ...
    'interfererCandidateFallback', interfererFallback, ...
    'peakSelection', peakDiagnostics, ...
    'usedTruth', false);
end

function [node, status] = strongestFiniteNode(powerMap, spectrumDiagnostics, mask)
%STRONGESTFINITENODE 为无峰或全零谱提供确定性的有限候选回退。
indices = find(mask & isfinite(powerMap));
if isempty(indices)
    node = zeros(0, 2);
    status = "unavailable";
    return;
end
[~, order] = max(real(powerMap(indices)));
linearIndex = indices(order);
node = [spectrumDiagnostics.xiX(linearIndex), ...
    spectrumDiagnostics.xiY(linearIndex)];
status = "available";
end

function [nodes, wasAdded] = appendUniqueNode(nodes, node, tolerance)
%APPENDUNIQUENODE 按配置容差保留首次候选，避免字典出现重复列。
wasAdded = false;
if isempty(node)
    return;
end
if all(max(abs(nodes - node), [], 2) > tolerance)
    nodes(end + 1, :) = node;
    wasAdded = true;
end
end

function nodes = uniqueNodes(nodes, existingNodes, tolerance)
%UNIQUENODES 保持峰值顺序并移除与已有候选近重复的方向节点。
uniqueList = zeros(0, 2);
for nodeIndex = 1:size(nodes, 1)
    node = nodes(nodeIndex, :);
    if isempty(existingNodes) || ...
            all(max(abs(existingNodes - node), [], 2) > tolerance)
        [uniqueList, wasAdded] = appendUniqueNode(uniqueList, node, tolerance);
        if wasAdded
            existingNodes = [existingNodes; node]; %#ok<AGROW>
        end
    end
end
nodes = uniqueList;
end

function [noisePower, source, diagnostics] = estimateNoisePower(Rx, snapshotCount, cfg)
%ESTIMATENOISEPOWER 按实际快拍数选择配置回退或小特征值稳健估计。
M = cfg.derived.M;
if snapshotCount < M
    noisePower = cfg.signal.noisePower;
    source = "configured-underdetermined";
    diagnostics = struct('eigenvalues', zeros(0, 1), ...
        'lowerEigenvalueCount', 0);
    return;
end
eigenvalues = sort(real(eig((Rx + Rx') / 2)), 'ascend');
lowerCount = max(1, ceil(M * cfg.spectrum.noiseEigenvalueFraction));
rawEstimate = median(eigenvalues(1:lowerCount));
scale = max(real(trace(Rx)) / M, eps);
noisePower = max(rawEstimate, eps(scale));
source = "lower-eigenvalues";
diagnostics = struct('eigenvalues', eigenvalues, ...
    'lowerEigenvalueCount', lowerCount);
end

function validateInputs(obs, cfg)
%VALIDATEINPUTS 校验白名单观测和CovFit配置，不接收任何真值字段。
if ~isstruct(obs) || isfield(obs, 'truth') || ~isfield(obs, 'X') || ...
        ~isfield(obs, 'Rx') || ~ismatrix(obs.X) || isempty(obs.X) || ...
        any(~isfinite(obs.X), 'all')
    error('uav:incm:InvalidObservation', ...
        'obs must be a truth-free observation with finite snapshots and covariance.');
end
if ~isstruct(cfg) || ~isfield(cfg, 'derived') || ...
        ~isfield(cfg.derived, 'M') || ...
        ~isfield(cfg.derived, 'nominalDirectionCosines') || ...
        ~isfield(cfg, 'incm') || ~isfield(cfg.incm, 'staticCovfit') || ...
        ~isfield(cfg, 'signal') || ~isfield(cfg, 'spectrum') || ...
        ~isfield(cfg, 'array')
    error('uav:incm:InvalidConfiguration', ...
        'cfg is missing UPA covariance-fitting fields.');
end
M = cfg.derived.M;
if size(obs.X, 1) ~= M || ~isequal(size(obs.Rx), [M, M]) || ...
        any(~isfinite(obs.Rx), 'all')
    error('uav:incm:InvalidObservation', ...
        'obs snapshots and covariance must match the configured array size.');
end
end
