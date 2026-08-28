function [RinHat, nextState, diagnostics] = reconstruct(method, obs, priorState, cfg)
%RECONSTRUCT 静态及时序INCM重构的统一入口。
%   method     - "static"或"temporal"；
%   obs        - 非Oracle白名单观测，至少含M×M样本协方差Rx；
%   priorState - 上一周期的估计谱与对齐诊断，静态方法可传[]；
%   RinHat     - M×M重构干扰加噪声协方差；
%   nextState  - 允许传递到下一周期的算法状态；
%   diagnostics- 谱、噪声估计、求积和结构投影等诊断信息。
%
%   本文件不访问真实方向、真实INCM或真实干扰数J。

method = lower(string(method));
validateInputs(obs, cfg);
switch method
    case "static"
        [RinHat, nextState, diagnostics] = reconstructStatic(obs, cfg);
    case "temporal"
        [RinHat, nextState, diagnostics] = reconstructTemporal(obs, priorState, cfg);
    otherwise
        error('uav:incm:UnsupportedMethod', ...
            'Supported INCM reconstruction methods are "static" and "temporal".');
end
end

function [RinHat, nextState, diagnostics] = reconstructStatic(obs, cfg)
% 静态基线：当前周期Capon谱 -> 期望区排除 -> 均匀求积 -> 结构投影。
[prepared, commonDiagnostics] = prepareSpectrum(obs, cfg, false);
quadrature = uniformQuadrature(prepared.residualPowerMap, ...
    prepared.spectrum, prepared.includedMask);
rawIncm = reconstructFromQuadrature(quadrature, prepared.noiseEstimate, cfg);
[RinHat, structureDiagnostics] = uav.incm.projectStructure(rawIncm, cfg, true);
nextState = struct('method', "static");
diagnostics = commonOutputDiagnostics(prepared, commonDiagnostics, ...
    structureDiagnostics);
diagnostics.usedPowerMap = prepared.residualPowerMap;
diagnostics.eigenvalueFloor = structureDiagnostics.eigenvalueFloor;
end

function [RinHat, nextState, diagnostics] = reconstructTemporal(obs, priorState, cfg)
% 步骤1：可选训练序列投影，并得到去噪声底座后的当前候选谱。
applyPilotProjection = cfg.incm.ablation.pilotProjection;
[prepared, commonDiagnostics] = prepareSpectrum(obs, cfg, applyPilotProjection);
% 步骤2：数据驱动估计上一周期谱到当前谱的二维平移。
[alignedPrevious, alignmentDiagnostics] = alignPrior( ...
    prepared.residualPowerMap, prepared.includedMask, priorState, cfg);
if cfg.incm.ablation.temporalContinuity
    effectiveLambdaT = cfg.incm.lambdaT;
else
    effectiveLambdaT = 0;
end
% 步骤3：仅在允许的干扰角域内求非负时序稀疏近端解。
enhancedSpectrum = zeros(size(prepared.residualPowerMap));
enhancedSpectrum(prepared.includedMask) = uav.incm.temporalProx( ...
    prepared.residualPowerMap(prepared.includedMask), ...
    alignedPrevious(prepared.includedMask), cfg.incm.lambdaS, effectiveLambdaT);
[peaks, peakDiagnostics] = uav.spectrum.selectPeaks(enhancedSpectrum, ...
    prepared.spectrum, prepared.includedMask, cfg);
% 步骤4：主算法使用峰值邻域细、其余区域粗的自适应求积；消融使用均匀网格。
if cfg.incm.ablation.adaptiveQuadrature
    quadrature = uav.spectrum.adaptiveQuadrature(enhancedSpectrum, ...
        prepared.spectrum, prepared.includedMask, peakDiagnostics.linearIndices, cfg);
    quadratureMode = "adaptive";
else
    quadrature = uniformQuadrature(enhancedSpectrum, ...
        prepared.spectrum, prepared.includedMask);
    quadratureMode = "uniform";
end
% 步骤5：积分得到原始INCM，再执行可消融的结构约束。
rawIncm = reconstructFromQuadrature(quadrature, prepared.noiseEstimate, cfg);
[RinHat, structureDiagnostics] = uav.incm.projectStructure(rawIncm, cfg, ...
    cfg.incm.ablation.structureEnhancement);

fallbacks = struct( ...
    'alignment', alignmentDiagnostics.fallbackReason, ...
    'noiseEstimate', commonDiagnostics.noise.source, ...
    'quadrature', quadratureMode);
% nextState中没有任何真实方向或真实协方差。
nextState = struct( ...
    'enhancedSpectrum', enhancedSpectrum, ...
    'xiX', prepared.spectrum.xiX, ...
    'xiY', prepared.spectrum.xiY, ...
    'includedMask', prepared.includedMask, ...
    'alignment', alignmentDiagnostics, ...
    'diagnostics', fallbacks);
diagnostics = commonOutputDiagnostics(prepared, commonDiagnostics, ...
    structureDiagnostics);
diagnostics.usedPowerMap = enhancedSpectrum;
diagnostics.enhancedSpectrum = enhancedSpectrum;
diagnostics.alignedPreviousSpectrum = alignedPrevious;
diagnostics.effectiveLambdaT = effectiveLambdaT;
diagnostics.peaks = peaks;
diagnostics.peakSelection = peakDiagnostics;
diagnostics.nodes = quadrature.nodes;
diagnostics.weights = quadrature.weights;
diagnostics.quadrature = quadrature.diagnostics;
diagnostics.quadratureMode = quadratureMode;
diagnostics.alignment = alignmentDiagnostics;
diagnostics.fallbacks = fallbacks;
diagnostics.rawIncm = rawIncm;
diagnostics.eigenvalueFloor = structureDiagnostics.eigenvalueFloor;
end

function [prepared, diagnostics] = prepareSpectrum(obs, cfg, applyPilotProjection)
% 构造两种INCM方法共享的Capon谱、允许掩码、噪声估计和残差谱。
[spectrumInputCovariance, projectionDiagnostics] = spectrumCovariance( ...
    obs, cfg, applyPilotProjection);
spectrumObservation = obs;
spectrumObservation.Rx = spectrumInputCovariance;
[powerMap, spectrumDiagnostics] = uav.spectrum.loadedCapon( ...
    spectrumInputCovariance, cfg);
includedMask = nonSoiVisibleMask(spectrumDiagnostics, cfg);
[noiseEstimate, noiseDiagnostics] = estimateNoisePower(spectrumObservation, cfg);
M = cfg.derived.M;
% 纯白噪声在单位幅度导向矢量下产生近似(noise+loading)/M的Capon底座。
pedestalPower = (noiseEstimate + spectrumDiagnostics.loading) / M;
residualPowerMap = zeros(size(powerMap));
residualPowerMap(spectrumDiagnostics.visibleMask) = max( ...
    powerMap(spectrumDiagnostics.visibleMask) - pedestalPower, 0);
prepared = struct( ...
    'powerMap', powerMap, ...
    'residualPowerMap', residualPowerMap, ...
    'spectrum', spectrumDiagnostics, ...
    'includedMask', includedMask, ...
    'noiseEstimate', noiseEstimate, ...
    'pedestalPower', pedestalPower, ...
    'spectrumInputCovariance', spectrumInputCovariance);
diagnostics = struct('noise', noiseDiagnostics, 'projection', projectionDiagnostics);
end

function [Rx, diagnostics] = spectrumCovariance(obs, cfg, applyProjection)
if ~applyProjection
    Rx = obs.Rx;
    diagnostics = struct('applied', false, 'removedPilotEnergy', 0);
    return;
end
if ~isfield(obs, 'X') || ~isfield(obs, 'pilot') || isempty(obs.X) || isempty(obs.pilot) || ...
        size(obs.X, 2) ~= numel(obs.pilot)
    error('uav:incm:MissingPilotData', ...
        'Pilot projection requires matching obs.X and obs.pilot samples.');
end
pilot = reshape(obs.pilot, 1, []);
pilotEnergy = real(pilot * pilot');
if ~isfinite(pilotEnergy) || pilotEnergy <= 0
    error('uav:incm:InvalidPilot', 'obs.pilot must have positive finite energy.');
end
% 将每个阵元快拍投影到已知训练序列张成的一维时间子空间。
pilotComponent = (obs.X * pilot') / pilotEnergy * pilot;
projectedSnapshots = obs.X - pilotComponent;
Rx = projectedSnapshots * projectedSnapshots' / size(projectedSnapshots, 2);
Rx = (Rx + Rx') / 2;
diagnostics = struct( ...
    'applied', true, ...
    'removedPilotEnergy', norm(pilotComponent, 'fro')^2 / size(projectedSnapshots, 2), ...
    'snapshotCount', size(projectedSnapshots, 2), ...
    'arraySize', cfg.derived.M);
end

function [alignedPrevious, diagnostics] = alignPrior(currentMap, includedMask, priorState, cfg)
% 首周期无历史状态时以全零历史谱开始，不虚构时间先验。
if isempty(priorState) || ~isstruct(priorState) || ...
        ~isfield(priorState, 'enhancedSpectrum')
    alignedPrevious = zeros(size(currentMap));
    diagnostics = identityAlignment("no-prior");
    return;
end
previous = priorState.enhancedSpectrum;
if ~isequal(size(previous), size(currentMap)) || ...
        any(~isfinite(previous), 'all') || any(previous < 0, 'all')
    alignedPrevious = zeros(size(currentMap));
    diagnostics = identityAlignment("invalid-prior");
    return;
end
previous(~includedMask) = 0;
if ~cfg.incm.ablation.alignment
    % G_k=I消融：保留历史谱，但不估计平移。
    alignedPrevious = previous;
    diagnostics = identityAlignment("disabled");
    return;
end
[alignedPrevious, diagnostics] = uav.spectrum.align2D( ...
    currentMap, previous, includedMask, cfg);
end

function diagnostics = identityAlignment(reason)
diagnostics = struct( ...
    'shiftBins', [0, 0], ...
    'bestCorrelation', NaN, ...
    'secondCorrelation', NaN, ...
    'ambiguityGap', NaN, ...
    'usedIdentity', true, ...
    'fallbackReason', string(reason), ...
    'candidateShifts', zeros(0, 2), ...
    'candidateCorrelations', zeros(0, 1));
end

function quadrature = uniformQuadrature(powerMap, spectrumDiagnostics, includedMask)
% 每个允许网格单元使用相同面积权重，作为静态基线和求积消融。
linearIndices = find(includedMask);
quadrature = struct( ...
    'nodes', [spectrumDiagnostics.xiX(linearIndices), ...
        spectrumDiagnostics.xiY(linearIndices)], ...
    'linearIndices', linearIndices, ...
    'weights', spectrumDiagnostics.cellArea * ones(numel(linearIndices), 1), ...
    'values', powerMap(linearIndices), ...
    'diagnostics', struct( ...
        'mode', "uniform", ...
        'fineCellCount', 0, ...
        'coarseCellCount', numel(linearIndices), ...
        'fineNodeCount', 0, ...
        'coarseNodeCount', numel(linearIndices), ...
        'representedArea', numel(linearIndices) * spectrumDiagnostics.cellArea));
end

function R = reconstructFromQuadrature(quadrature, noiseEstimate, cfg)
% 实现sum_n weight_n*power_n*a_n*a_n^H + sigma_n^2*I。
M = cfg.derived.M;
if isfield(quadrature, 'linearIndices')
    % 均匀网格节点直接复用config阶段缓存的流形列。
    grid = cfg.derived.spectrumGrid;
    [found, cachedColumns] = ismember(quadrature.linearIndices, ...
        grid.visibleLinearIndices);
    if ~all(found)
        error('uav:incm:InvalidQuadrature', ...
            'Uniform quadrature nodes must belong to the cached visible grid.');
    end
    manifold = grid.manifold(:, cachedColumns);
else
    % 自适应细节点不在原网格上，但仍一次批量生成全部导向矢量。
    manifold = uav.array.steeringUPA(quadrature.nodes.', ...
        cfg.array.Nx, cfg.array.Ny, cfg.array.spacingOverLambda);
end
weightedPower = reshape(quadrature.weights .* quadrature.values, 1, []);
% (A.*p)*A'等价于A*diag(p)*A'，但不显式构造大型对角矩阵。
signalCovariance = (manifold .* weightedPower) * manifold';
R = signalCovariance + noiseEstimate * eye(M);
R = (R + R') / 2;
end

function diagnostics = commonOutputDiagnostics(prepared, commonDiagnostics, structureDiagnostics)
spectrumDiagnostics = prepared.spectrum;
noiseDiagnostics = commonDiagnostics.noise;
diagnostics = struct( ...
    'rawPowerMap', prepared.powerMap, ...
    'residualPowerMap', prepared.residualPowerMap, ...
    'powerMap', prepared.powerMap, ...
    'spectrum', spectrumDiagnostics, ...
    'includedMask', prepared.includedMask, ...
    'excludedSoiMask', spectrumDiagnostics.visibleMask & ~prepared.includedMask, ...
    'noiseEstimate', prepared.noiseEstimate, ...
    'noiseFloor', noiseDiagnostics.noiseFloor, ...
    'effectiveSnapshotCount', noiseDiagnostics.effectiveSnapshotCount, ...
    'effectiveRank', noiseDiagnostics.effectiveRank, ...
    'noiseEstimateSource', noiseDiagnostics.source, ...
    'pedestalPower', prepared.pedestalPower, ...
    'structure', structureDiagnostics, ...
    'pilotProjectionApplied', commonDiagnostics.projection.applied, ...
    'pilotProjection', commonDiagnostics.projection, ...
    'spectrumInputCovariance', prepared.spectrumInputCovariance);
end

function validateInputs(obs, cfg)
if ~isstruct(obs) || ~isfield(obs, 'Rx')
    error('uav:incm:InvalidObservation', 'obs must provide the sample covariance Rx.');
end
if ~isstruct(cfg) || ~isfield(cfg, 'derived') || ~isfield(cfg, 'incm') || ...
        ~isfield(cfg, 'spectrum') || ~isfield(cfg, 'array')
    error('uav:incm:InvalidConfiguration', 'cfg is missing INCM fields.');
end
M = cfg.derived.M;
if ~isequal(size(obs.Rx), [M, M]) || any(~isfinite(obs.Rx), 'all')
    error('uav:incm:InvalidObservation', ...
        'obs.Rx must be finite and match the configured array size.');
end
if isfield(obs, 'X') && (~ismatrix(obs.X) || size(obs.X, 1) ~= M || ...
        size(obs.X, 2) < 1 || any(~isfinite(obs.X), 'all'))
    error('uav:incm:InvalidObservation', ...
        'obs.X must be finite, nonempty, and have the configured sensor rows.');
end
end

function includedMask = nonSoiVisibleMask(spectrumDiagnostics, cfg)
% 同时施加可见圆盘与连续三维期望信号排除区约束。
xiX = spectrumDiagnostics.xiX;
xiY = spectrumDiagnostics.xiY;
continuousMask = uav.spectrum.admissibleDirections([xiX(:), xiY(:)], cfg);
includedMask = spectrumDiagnostics.visibleMask & reshape(continuousMask, size(xiX));
end

function [noiseEstimate, diagnostics] = estimateNoisePower(obs, cfg)
% 使用SCM较小非零特征值的稳健统计估计噪声，不读取真实干扰数。
Rx = obs.Rx;
eigenvalues = sort(real(eig((Rx + Rx') / 2)), 'ascend');
M = numel(eigenvalues);
covarianceScale = max(real(trace(Rx)) / M, eps);
noiseFloor = cfg.incm.noiseFloorRatio * covarianceScale;
[snapshotCount, hasSnapshots] = observedSnapshotCount(obs);
zeroTolerance = eps(max(1, max(abs(eigenvalues)))) * max(M, snapshotCount);
nonzeroEigenvalues = eigenvalues(eigenvalues > zeroTolerance);
if ~hasSnapshots
    % 只有协方差而无原始快拍时，无法判断有效快拍秩，使用配置接收机噪声。
    noiseEstimate = cfg.signal.noisePower;
    source = "configured-missing-snapshots";
elseif snapshotCount < M
    % L<M时SCM必然秩亏，零特征值不能作为噪声功率估计，按合同回退。
    noiseEstimate = cfg.signal.noisePower;
    source = "configured-underdetermined";
elseif isempty(nonzeroEigenvalues)
    noiseEstimate = noiseFloor;
    source = "lower-eigenvalues";
else
    % L>=M时取较小非零特征值中的配置比例，并用中位数抑制干扰离群值。
    lowerCount = max(1, ceil(cfg.spectrum.noiseEigenvalueFraction * ...
        numel(nonzeroEigenvalues)));
    lowerEstimate = median(nonzeroEigenvalues(1:lowerCount));
    noiseEstimate = max(lowerEstimate, noiseFloor);
    source = "lower-eigenvalues";
end
diagnostics = struct( ...
    'noiseFloor', noiseFloor, ...
    'effectiveSnapshotCount', snapshotCount, ...
    'effectiveRank', numel(nonzeroEigenvalues), ...
    'source', source);
end

function [snapshotCount, hasSnapshots] = observedSnapshotCount(obs)
hasSnapshots = isfield(obs, 'X') && ~isempty(obs.X);
if hasSnapshots
    snapshotCount = size(obs.X, 2);
else
    snapshotCount = 0;
end
end
