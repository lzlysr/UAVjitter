function [RinHat, nextState, diagnostics] = reconstruct(method, obs, priorState, cfg)
%RECONSTRUCT Estimate a static interference-plus-noise covariance matrix.

method = lower(string(method));
if method ~= "static"
    error('uav:incm:UnsupportedMethod', ...
        'Task 3 supports only the "static" INCM reconstruction method.');
end
validateInputs(obs, cfg);
priorState = priorState; %#ok<NASGU>

[powerMap, spectrumDiagnostics] = uav.spectrum.loadedCapon(obs.Rx, cfg);
includedMask = nonSoiVisibleMask(spectrumDiagnostics, cfg);
M = cfg.derived.M;
signalCovariance = zeros(M);
for linearIndex = find(includedMask).'
    a = uav.array.steeringUPA([spectrumDiagnostics.xiX(linearIndex); ...
        spectrumDiagnostics.xiY(linearIndex)], cfg.array.Nx, cfg.array.Ny, ...
        cfg.array.spacingOverLambda);
    signalCovariance = signalCovariance + spectrumDiagnostics.cellArea * ...
        powerMap(linearIndex) * (a * a');
end

noiseEstimate = estimateNoisePower(obs.Rx, cfg);
rawIncm = signalCovariance + noiseEstimate * eye(M);
[RinHat, eigenvalueFloor] = enforceStructure(rawIncm, cfg);
nextState = struct('method', "static");
diagnostics = struct( ...
    'powerMap', powerMap, ...
    'spectrum', spectrumDiagnostics, ...
    'includedMask', includedMask, ...
    'excludedSoiMask', spectrumDiagnostics.visibleMask & ~includedMask, ...
    'noiseEstimate', noiseEstimate, ...
    'eigenvalueFloor', eigenvalueFloor);
end

function validateInputs(obs, cfg)
if ~isstruct(obs) || ~isfield(obs, 'Rx')
    error('uav:incm:InvalidObservation', 'obs must provide the sample covariance Rx.');
end
if ~isstruct(cfg) || ~isfield(cfg, 'derived') || ~isfield(cfg, 'incm') || ...
        ~isfield(cfg, 'spectrum') || ~isfield(cfg, 'array')
    error('uav:incm:InvalidConfiguration', 'cfg is missing static INCM fields.');
end
M = cfg.derived.M;
if ~isequal(size(obs.Rx), [M, M]) || any(~isfinite(obs.Rx), 'all')
    error('uav:incm:InvalidObservation', ...
        'obs.Rx must be finite and match the configured array size.');
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

function noiseEstimate = estimateNoisePower(Rx, cfg)
eigenvalues = sort(real(eig((Rx + Rx') / 2)), 'ascend');
M = numel(eigenvalues);
lowerCount = max(1, ceil(cfg.spectrum.noiseEigenvalueFraction * M));
lowerEstimate = median(eigenvalues(1:lowerCount));
covarianceScale = max(real(trace(Rx)) / M, eps);
noiseEstimate = max(lowerEstimate, cfg.incm.noiseFloorRatio * covarianceScale);
end

function [structured, eigenvalueFloor] = enforceStructure(R, cfg)
M = size(R, 1);
R = (R + R') / 2;
[vectors, values] = eig(R, 'vector');
values = real(values);
covarianceScale = max(real(trace(R)) / M, eps);
eigenvalueFloor = cfg.incm.eigenvalueFloorRatio * covarianceScale;
values = max(values, eigenvalueFloor);
structured = vectors * diag(values) * vectors';
structured = (structured + structured') / 2;
isotropicCovariance = real(trace(structured)) / M * eye(M);
structured = (1 - cfg.incm.shrinkage) * structured + ...
    cfg.incm.shrinkage * isotropicCovariance;
structured = (structured + structured') / 2;
end
