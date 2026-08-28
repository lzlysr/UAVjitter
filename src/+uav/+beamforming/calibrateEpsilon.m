function [epsilonS, diagnostics] = calibrateEpsilon(cfg, stream)
%CALIBRATEEPSILON Calibrate a fixed steering-error radius from design jitter.

if nargin < 2 || isempty(stream)
    stream = uav.util.makeStream(cfg.random.baseSeed, ...
        cfg.beamforming.calibration.substream);
end
validateInputs(cfg, stream);
calibration = cfg.beamforming.calibration;
sampleCount = calibration.sampleCount;
designStdRad = deg2rad(reshape(calibration.designJitterStdDeg, 3, 1));
offsets = designStdRad .* randn(stream, 3, sampleCount);
nominalAttitudeRad = deg2rad(reshape(cfg.attitude.nominalDeg, 3, 1));
desiredWorldDirection = cfg.derived.worldUnitVectors(:, 1);
aRef = cfg.derived.nominalSteeringReference;
mismatchNorms = zeros(sampleCount, 1);
for index = 1:sampleCount
    rotation = uav.geometry.rotationZYX(nominalAttitudeRad + offsets(:, index));
    arrayDirection = uav.geometry.worldToArray(rotation, desiredWorldDirection);
    directionCosines = uav.geometry.directionCoordinates(arrayDirection);
    steering = uav.array.steeringUPA(directionCosines, cfg.array.Nx, ...
        cfg.array.Ny, cfg.array.spacingOverLambda);
    mismatchNorms(index) = norm(steering - aRef);
end
epsilonS = quantile(mismatchNorms, calibration.quantile);
if ~isfinite(epsilonS) || epsilonS < 0 || epsilonS >= norm(aRef)
    error('uav:beamforming:InfeasibleCalibration', ...
        'The calibrated epsilon must be finite and smaller than ||aRef||.');
end
diagnostics = struct( ...
    'designJitterStdDeg', calibration.designJitterStdDeg, ...
    'quantile', calibration.quantile, ...
    'sampleCount', sampleCount, ...
    'configuredSubstream', calibration.substream, ...
    'mismatchNorms', mismatchNorms, ...
    'referenceNorm', norm(aRef));
end

function validateInputs(cfg, stream)
if ~isstruct(cfg) || ~isfield(cfg, 'beamforming') || ...
        ~isfield(cfg.beamforming, 'calibration') || ~isfield(cfg, 'derived') || ...
        ~isfield(cfg, 'attitude') || ~isfield(cfg, 'array')
    error('uav:beamforming:InvalidConfiguration', ...
        'cfg is missing offline calibration fields.');
end
if ~isa(stream, 'RandStream')
    error('uav:beamforming:InvalidStream', ...
        'Calibration requires an explicit RandStream.');
end
end
