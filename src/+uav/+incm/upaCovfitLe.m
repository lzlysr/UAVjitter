function [RinHat, diagnostics] = upaCovfitLe(obs, cfg)
%   移植自Gu–Leshem研究分支d666399；仅适配配置路径及函数名。
%UPACOVFIT 以加载Capon候选和非负协方差拟合重构UPA干扰协方差。
%   该方法只使用算法白名单中的快拍、样本协方差和配置先验；真实方向、
%   真实姿态、真实干扰数及真实INCM均不参与候选、噪声或功率计算。

validateInputs(obs, cfg);
M = cfg.derived.M;
diagnosis = cfg.experiments.linYueStudy.leCovfit;
snapshotCount = size(obs.X, 2);

% 仅在局部副本中映射诊断方法的加载和峰值参数，原配置保持不变。
spectrumCfg = cfg;
spectrumCfg.spectrum.diagonalLoading = diagnosis.loadingRatio;
spectrumCfg.spectrum.peakThresholdRatio = diagnosis.peakThresholdRatio;
spectrumCfg.spectrum.peakNmsRadiusBins = diagnosis.peakNmsRadiusBins;
[powerMap, spectrumDiagnostics] = uav.spectrum.loadedCapon(obs.Rx, spectrumCfg);

allNodes = [spectrumDiagnostics.xiX(:), spectrumDiagnostics.xiY(:)];
outsideDesired = uav.spectrum.admissibleDirections(allNodes, cfg);
interfererMask = spectrumDiagnostics.visibleMask & ...
    reshape(outsideDesired, size(spectrumDiagnostics.visibleMask));
desiredMask = spectrumDiagnostics.visibleMask & ~interfererMask;
nominalNode = cfg.derived.nominalDirectionCosines(:, 1).';
tolerance = cfg.spectrum.quadratureDedupTolerance;

% 期望字典始终含名义方向；期望区最强有限谱点只在非重复时加入。
[desiredPeakNode, desiredPeakStatus] = strongestFiniteNode( ...
    powerMap, spectrumDiagnostics, desiredMask);
[desiredNodes, desiredPeakAdded] = appendUniqueNode(nominalNode, ...
    desiredPeakNode, tolerance);
if desiredPeakStatus == "unavailable"
    desiredFallback = "nominal-only";
elseif desiredPeakAdded
    desiredFallback = "strongest-finite";
else
    desiredFallback = "deduplicated-with-nominal";
end

% 干扰区复用既有阈值与二维NMS语义，再按本方法的候选数上限截断。
[peakNodes, peakDiagnostics] = uav.spectrum.selectPeaks( ...
    powerMap, spectrumDiagnostics, interfererMask, spectrumCfg);
candidateCap = min(diagnosis.maximumCandidateCount, M - 1);
peakNodes = peakNodes(1:min(size(peakNodes, 1), candidateCap), :);
interfererFallback = "nms";
if candidateCap == 0
    % 单阵元没有可容纳的干扰空间模态，不能由无峰回退突破硬候选上限。
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

% 专项根因研究的显式候选细化开关；默认关闭，既有方法定义保持不变。
% 开启后必须另列方法名称：这不是原始Rec，也不是历史CovFit结果。
refinementDiagnostics = struct('enabled', false);
if isfield(diagnosis, 'refineCandidates') && diagnosis.refineCandidates
    [desiredRefined, desiredRefinement] = uav.spectrum.refineCaponNodesLe( ...
        desiredNodes(2:end, :), spectrumDiagnostics.loadedCovariance, true, cfg);
    [interfererRefined, interfererRefinement] = uav.spectrum.refineCaponNodesLe( ...
        interfererNodes, spectrumDiagnostics.loadedCovariance, false, cfg);
    refinedTolerance = cfg.experiments.linYueStudy.leRefinement.dedupTolerance;
    desiredNodes = uniqueNodes([nominalNode; desiredRefined], zeros(0, 2), refinedTolerance);
    interfererNodes = uniqueNodes(interfererRefined, desiredNodes, refinedTolerance);
    refinementDiagnostics = struct('enabled', true, ...
        'desired', desiredRefinement, 'interferer', interfererRefinement);
end

[noisePower, noiseSource, noiseDiagnostics] = estimateNoisePower( ...
    obs.Rx, snapshotCount, cfg);
candidateNodes = [desiredNodes; interfererNodes];
candidateManifold = uav.array.steeringUPA(candidateNodes.', cfg.array.Nx, ...
    cfg.array.Ny, cfg.array.spacingOverLambda);
[fittedPowers, fitDiagnostics] = uav.incm.fitNonnegativePowers( ...
    obs.Rx, candidateManifold, noisePower, diagnosis.fitRegularization);

% 最终INCM仅保留干扰候选与独立噪声项，期望候选只服务于联合拟合。
interfererPowers = fittedPowers(size(desiredNodes, 1) + 1:end);
interfererManifold = candidateManifold(:, size(desiredNodes, 1) + 1:end);
if isempty(interfererPowers)
    interferenceCovariance = zeros(M);
else
    interferenceCovariance = (interfererManifold .* interfererPowers.') * ...
        interfererManifold';
end
RinHat = interferenceCovariance + noisePower * eye(M);
RinHat = (RinHat + RinHat') / 2;

diagnostics = struct( ...
    'candidateRefinement', refinementDiagnostics, ...
    'loading', spectrumDiagnostics.loading, ...
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
    'desiredFittedPowers', fittedPowers(1:size(desiredNodes, 1)), ...
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
values = real(powerMap(indices));
[~, order] = max(values);
linearIndex = indices(order);
node = [spectrumDiagnostics.xiX(linearIndex), ...
    spectrumDiagnostics.xiY(linearIndex)];
status = "available";
end

function [nodes, wasAdded] = appendUniqueNode(nodes, node, tolerance)
%APPENDUNIQUENODE 按配置容差保留首次候选，避免字典重复列。
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
%UNIQUENODES 保持NMS顺序并移除与已有候选近重复的方向节点。
uniqueNodes = zeros(0, 2);
for nodeIndex = 1:size(nodes, 1)
    node = nodes(nodeIndex, :);
    if isempty(existingNodes) || ...
            all(max(abs(existingNodes - node), [], 2) > tolerance)
        [uniqueNodes, wasAdded] = appendUniqueNode(uniqueNodes, node, tolerance);
        if wasAdded
            existingNodes = [existingNodes; node]; %#ok<AGROW>
        end
    end
end
nodes = uniqueNodes;
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
%VALIDATEINPUTS 只校验算法白名单所需字段与协方差拟合配置。
if ~isstruct(obs) || isfield(obs, 'truth') || ~isfield(obs, 'X') || ~isfield(obs, 'Rx') || ...
        ~ismatrix(obs.X) || isempty(obs.X) || any(~isfinite(obs.X), 'all')
    error('uav:incm:InvalidObservation', ...
        'obs must provide finite nonempty snapshots X and covariance Rx.');
end
if ~isstruct(cfg) || ~isfield(cfg, 'derived') || ~isfield(cfg.derived, 'M') || ...
        ~isfield(cfg.derived, 'nominalDirectionCosines') || ...
        ~isfield(cfg, 'experiments') || ~isfield(cfg.experiments, 'linYueStudy') || ...
        ~isfield(cfg.experiments.linYueStudy, 'leCovfit') || ...
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
