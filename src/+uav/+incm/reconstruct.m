function [RinHat, nextState, diagnostics] = reconstruct(method, obs, priorState, cfg)
%RECONSTRUCT Estimate static or temporal interference-plus-noise covariance.

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
applyPilotProjection = cfg.incm.ablation.pilotProjection;
[prepared, commonDiagnostics] = prepareSpectrum(obs, cfg, applyPilotProjection);
[alignedPrevious, alignmentDiagnostics] = alignPrior( ...
    prepared.residualPowerMap, prepared.includedMask, priorState, cfg);
if cfg.incm.ablation.temporalContinuity
    effectiveLambdaT = cfg.incm.lambdaT;
else
    effectiveLambdaT = 0;
end
enhancedSpectrum = zeros(size(prepared.residualPowerMap));
enhancedSpectrum(prepared.includedMask) = uav.incm.temporalProx( ...
    prepared.residualPowerMap(prepared.includedMask), ...
    alignedPrevious(prepared.includedMask), cfg.incm.lambdaS, effectiveLambdaT);
[peaks, peakDiagnostics] = uav.spectrum.selectPeaks(enhancedSpectrum, ...
    prepared.spectrum, prepared.includedMask, cfg);
if cfg.incm.ablation.adaptiveQuadrature
    quadrature = uav.spectrum.adaptiveQuadrature(enhancedSpectrum, ...
        prepared.spectrum, prepared.includedMask, peakDiagnostics.linearIndices, cfg);
    quadratureMode = "adaptive";
else
    quadrature = uniformQuadrature(enhancedSpectrum, ...
        prepared.spectrum, prepared.includedMask);
    quadratureMode = "uniform";
end
rawIncm = reconstructFromQuadrature(quadrature, prepared.noiseEstimate, cfg);
[RinHat, structureDiagnostics] = uav.incm.projectStructure(rawIncm, cfg, ...
    cfg.incm.ablation.structureEnhancement);

fallbacks = struct( ...
    'alignment', alignmentDiagnostics.fallbackReason, ...
    'noiseEstimate', commonDiagnostics.noise.source, ...
    'quadrature', quadratureMode);
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
[spectrumInputCovariance, projectionDiagnostics] = spectrumCovariance( ...
    obs, cfg, applyPilotProjection);
spectrumObservation = obs;
spectrumObservation.Rx = spectrumInputCovariance;
[powerMap, spectrumDiagnostics] = uav.spectrum.loadedCapon( ...
    spectrumInputCovariance, cfg);
includedMask = nonSoiVisibleMask(spectrumDiagnostics, cfg);
[noiseEstimate, noiseDiagnostics] = estimateNoisePower(spectrumObservation, cfg);
M = cfg.derived.M;
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
    'usedIdentity', true, ...
    'fallbackReason', string(reason), ...
    'candidateShifts', zeros(0, 2), ...
    'candidateCorrelations', zeros(0, 1));
end

function quadrature = uniformQuadrature(powerMap, spectrumDiagnostics, includedMask)
linearIndices = find(includedMask);
quadrature = struct( ...
    'nodes', [spectrumDiagnostics.xiX(linearIndices), ...
        spectrumDiagnostics.xiY(linearIndices)], ...
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
M = cfg.derived.M;
signalCovariance = zeros(M);
for index = 1:size(quadrature.nodes, 1)
    a = uav.array.steeringUPA(quadrature.nodes(index, :).', ...
        cfg.array.Nx, cfg.array.Ny, cfg.array.spacingOverLambda);
    signalCovariance = signalCovariance + quadrature.weights(index) * ...
        quadrature.values(index) * (a * a');
end
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
xiX = spectrumDiagnostics.xiX;
xiY = spectrumDiagnostics.xiY;
z = -sqrt(max(0, 1 - xiX.^2 - xiY.^2));
nominalDirection = cfg.derived.nominalArrayUnitVectors(:, 1);
cosine = xiX * nominalDirection(1) + xiY * nominalDirection(2) + ...
    z * nominalDirection(3);
angularDistance = acos(max(-1, min(1, real(cosine))));
isSoi = angularDistance <= deg2rad(cfg.spectrum.desiredExclusionHalfWidthDeg);
includedMask = spectrumDiagnostics.visibleMask & ~isSoi;
end

function [noiseEstimate, diagnostics] = estimateNoisePower(obs, cfg)
Rx = obs.Rx;
eigenvalues = sort(real(eig((Rx + Rx') / 2)), 'ascend');
M = numel(eigenvalues);
covarianceScale = max(real(trace(Rx)) / M, eps);
noiseFloor = cfg.incm.noiseFloorRatio * covarianceScale;
[snapshotCount, hasSnapshots] = observedSnapshotCount(obs);
zeroTolerance = eps(max(1, max(abs(eigenvalues)))) * max(M, snapshotCount);
nonzeroEigenvalues = eigenvalues(eigenvalues > zeroTolerance);
if ~hasSnapshots
    noiseEstimate = cfg.signal.noisePower;
    source = "configured-missing-snapshots";
elseif snapshotCount < M
    noiseEstimate = cfg.signal.noisePower;
    source = "configured-underdetermined";
elseif isempty(nonzeroEigenvalues)
    noiseEstimate = noiseFloor;
    source = "lower-eigenvalues";
else
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
