function cfg = validate(cfg)
%VALIDATE Check simulator configuration and derive array/scenario quantities.

requiredSections = {'random', 'array', 'timeline', 'attitude', 'sources', ...
    'signal', 'spectrum', 'incm', 'beamforming', 'monteCarlo', 'experiments', ...
    'plot', 'output'};
for index = 1:numel(requiredSections)
    if ~isfield(cfg, requiredSections{index})
        error('uav:config:MissingSection', 'Missing configuration section "%s".', ...
            requiredSections{index});
    end
end

validatePositiveInteger(cfg.array.Nx, 'array.Nx');
validatePositiveInteger(cfg.array.Ny, 'array.Ny');
validateFiniteVector(cfg.array.spacingOverLambda, 2, 'array.spacingOverLambda');
if any(cfg.array.spacingOverLambda <= 0)
    error('uav:config:NonPositiveSpacing', 'Array spacings must be positive.');
end
validatePositiveInteger(cfg.timeline.K, 'timeline.K');
validatePositiveInteger(cfg.timeline.L, 'timeline.L');
validateFiniteScalar(cfg.timeline.dtSec, 'timeline.dtSec');
if cfg.timeline.dtSec <= 0
    error('uav:config:NonPositiveDt', 'timeline.dtSec must be positive.');
end
validateFiniteVector(cfg.attitude.nominalDeg, 3, 'attitude.nominalDeg');
validateFiniteVector(cfg.attitude.stdDeg, 3, 'attitude.stdDeg');
if any(cfg.attitude.stdDeg < 0)
    error('uav:config:NegativeJitterStd', 'attitude.stdDeg must be nonnegative.');
end
validateFiniteVector(cfg.attitude.rho, 3, 'attitude.rho');
if any(abs(cfg.attitude.rho) >= 1)
    error('uav:config:InvalidRho', 'Each attitude correlation coefficient must lie in (-1, 1).');
end
validatePositiveInteger(cfg.random.baseSeed, 'random.baseSeed');
validatePositiveInteger(cfg.monteCarlo.nTrials, 'monteCarlo.nTrials');
validatePositiveInteger(cfg.spectrum.gridSize, 'spectrum.gridSize');
if cfg.spectrum.gridSize < 3 || mod(cfg.spectrum.gridSize, 2) == 0
    error('uav:config:InvalidPositiveInteger', ...
        'spectrum.gridSize must be an odd integer greater than or equal to three.');
end
validateAlgorithmParameters(cfg);
validateFiniteScalar(cfg.signal.snrDb, 'signal.snrDb');
validateFiniteScalar(cfg.signal.noisePower, 'signal.noisePower');
if cfg.signal.noisePower <= 0
    error('uav:config:NonPositiveNoisePower', 'signal.noisePower must be positive.');
end
validateFiniteVector(cfg.experiments.jitterSweepDeg, [], 'experiments.jitterSweepDeg');
validateFiniteVector(cfg.experiments.snapshotSweep, [], 'experiments.snapshotSweep');
if any(cfg.experiments.jitterSweepDeg < 0) || any(cfg.experiments.snapshotSweep <= 0)
    error('uav:config:InvalidSweep', 'Sweeps must contain nonnegative jitter and positive snapshots.');
end

validateDirections(cfg.sources.desiredAzElDeg, 'desired');
if ~isequal(size(cfg.sources.desiredAzElDeg), [1, 2])
    error('uav:config:DesiredDirectionShape', 'Desired direction must be a 1-by-2 azimuth/elevation vector.');
end
if isempty(cfg.sources.interfererAzElDeg)
    cfg.sources.interfererAzElDeg = zeros(0, 2);
end
if size(cfg.sources.interfererAzElDeg, 2) ~= 2
    error('uav:config:InterfererDirectionShape', ...
        'Interferer directions must have two columns: azimuth and elevation.');
end
validateDirections(cfg.sources.interfererAzElDeg, 'interferer');

J = size(cfg.sources.interfererAzElDeg, 1);
if J == 0 && isempty(cfg.sources.interfererInrDb)
    cfg.sources.interfererInrDb = zeros(1, 0);
else
    validateFiniteVector(cfg.sources.interfererInrDb, [], 'sources.interfererInrDb');
    cfg.sources.interfererInrDb = reshape(cfg.sources.interfererInrDb, 1, []);
end
if numel(cfg.sources.interfererInrDb) ~= J
    error('uav:config:InrLengthMismatch', ...
        'The number of interferer INR values must equal the number of interferer directions.');
end

azElDeg = [cfg.sources.desiredAzElDeg; cfg.sources.interfererAzElDeg];
worldUnitVectors = azElToUnitVectors(azElDeg);
nominalRotation = rotationZYX(deg2rad(cfg.attitude.nominalDeg));
nominalArrayUnitVectors = nominalRotation.' * worldUnitVectors;
directionCosines = nominalArrayUnitVectors(1:2, :);
nominalSteering = steeringReference(directionCosines, cfg.array.Nx, cfg.array.Ny, ...
    cfg.array.spacingOverLambda);

noisePower = cfg.signal.noisePower;
desiredPower = noisePower * 10^(cfg.signal.snrDb / 10);
interfererPowers = noisePower * 10.^(cfg.sources.interfererInrDb(:).' / 10);
rootDir = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
cfg.output = struct( ...
    'rootDir', rootDir, ...
    'cacheDir', fullfile(rootDir, 'data', 'cache'), ...
    'rawDir', fullfile(rootDir, 'results', 'raw'), ...
    'figureDir', fullfile(rootDir, 'results', 'figures'), ...
    'tableDir', fullfile(rootDir, 'results', 'tables'), ...
    'logDir', fullfile(rootDir, 'results', 'logs'));
cfg.derived = struct( ...
    'M', cfg.array.Nx * cfg.array.Ny, ...
    'J', J, ...
    'noisePower', noisePower, ...
    'desiredPower', desiredPower, ...
    'interfererPowers', interfererPowers, ...
    'sourcePowers', [desiredPower, interfererPowers], ...
    'worldUnitVectors', worldUnitVectors, ...
    'nominalRotation', nominalRotation, ...
    'nominalArrayUnitVectors', nominalArrayUnitVectors, ...
    'nominalDirectionCosines', directionCosines, ...
    'nominalSteering', nominalSteering, ...
    'nominalSteeringReference', nominalSteering(:, 1));
end

function validatePositiveInteger(value, name)
if ~isscalar(value) || ~isfinite(value) || value <= 0 || value ~= floor(value)
    error('uav:config:InvalidPositiveInteger', '%s must be a positive integer.', name);
end
end

function validateFiniteScalar(value, name)
if ~isscalar(value) || ~isfinite(value)
    error('uav:config:InvalidScalar', '%s must be a finite scalar.', name);
end
end

function validateFiniteVector(value, expectedLength, name)
if ~isvector(value) || ~isnumeric(value) || any(~isfinite(value)) || ...
        (~isempty(expectedLength) && numel(value) ~= expectedLength)
    error('uav:config:InvalidVector', '%s must be a finite vector of the required length.', name);
end
end

function validateAlgorithmParameters(cfg)
validateAlgorithmFiniteScalar(cfg.spectrum.desiredExclusionHalfWidthDeg, ...
    'spectrum.desiredExclusionHalfWidthDeg');
if cfg.spectrum.desiredExclusionHalfWidthDeg < 0 || ...
        cfg.spectrum.desiredExclusionHalfWidthDeg >= 90
    error('uav:config:InvalidAlgorithmParameter', ...
        'The desired exclusion half width must lie in [0, 90) degrees.');
end
validateNonnegativeFiniteScalar(cfg.spectrum.diagonalLoading, 'spectrum.diagonalLoading');
validatePositiveFiniteScalar(cfg.spectrum.noiseEigenvalueFraction, ...
    'spectrum.noiseEigenvalueFraction');
if cfg.spectrum.noiseEigenvalueFraction > 1
    error('uav:config:InvalidAlgorithmParameter', ...
        'spectrum.noiseEigenvalueFraction must not exceed one.');
end
validateAlgorithmFiniteScalar(cfg.incm.shrinkage, 'incm.shrinkage');
if cfg.incm.shrinkage < 0 || cfg.incm.shrinkage > 1
    error('uav:config:InvalidAlgorithmParameter', ...
        'incm.shrinkage must lie in [0, 1].');
end
validatePositiveFiniteScalar(cfg.incm.eigenvalueFloorRatio, ...
    'incm.eigenvalueFloorRatio');
if cfg.incm.eigenvalueFloorRatio > 1
    error('uav:config:InvalidAlgorithmParameter', ...
        'incm.eigenvalueFloorRatio must not exceed one.');
end
validatePositiveFiniteScalar(cfg.incm.noiseFloorRatio, 'incm.noiseFloorRatio');
if cfg.incm.noiseFloorRatio > 1
    error('uav:config:InvalidAlgorithmParameter', ...
        'incm.noiseFloorRatio must not exceed one.');
end
validateNonnegativeFiniteScalar(cfg.incm.lambdaS, 'incm.lambdaS');
validateNonnegativeFiniteScalar(cfg.incm.lambdaT, 'incm.lambdaT');
validateNonnegativeFiniteScalar(cfg.beamforming.diagonalLoading, ...
    'beamforming.diagonalLoading');
validateNonnegativeFiniteScalar(cfg.beamforming.epsilonS, 'beamforming.epsilonS');
validatePositiveFiniteScalar(cfg.beamforming.minimumReciprocalCondition, ...
    'beamforming.minimumReciprocalCondition');
if cfg.beamforming.minimumReciprocalCondition >= 1
    error('uav:config:InvalidAlgorithmParameter', ...
        'beamforming.minimumReciprocalCondition must lie in (0, 1).');
end
end

function validateNonnegativeFiniteScalar(value, name)
validateAlgorithmFiniteScalar(value, name);
if value < 0
    error('uav:config:InvalidAlgorithmParameter', '%s must be nonnegative.', name);
end
end

function validatePositiveFiniteScalar(value, name)
validateAlgorithmFiniteScalar(value, name);
if value <= 0
    error('uav:config:InvalidAlgorithmParameter', '%s must be positive.', name);
end
end

function validateAlgorithmFiniteScalar(value, name)
if ~isscalar(value) || ~isfinite(value)
    error('uav:config:InvalidAlgorithmParameter', ...
        '%s must be a finite scalar.', name);
end
end

function validateDirections(azElDeg, sourceName)
if ~isnumeric(azElDeg) || size(azElDeg, 2) ~= 2 || any(~isfinite(azElDeg), 'all')
    error('uav:config:InvalidDirection', '%s directions must be finite azimuth/elevation pairs.', sourceName);
end
elevationDeg = azElDeg(:, 2);
if any(elevationDeg < -90 | elevationDeg >= 0)
    error('uav:config:SourceNotLowerHemisphere', ...
        '%s elevations must be in the visible lower hemisphere [-90, 0).', sourceName);
end
end

function vectors = azElToUnitVectors(azElDeg)
azimuthRad = deg2rad(azElDeg(:, 1));
elevationRad = deg2rad(azElDeg(:, 2));
vectors = [cos(elevationRad) .* cos(azimuthRad), ...
    cos(elevationRad) .* sin(azimuthRad), sin(elevationRad)].';
end

function R = rotationZYX(etaRad)
roll = etaRad(1);
pitch = etaRad(2);
yaw = etaRad(3);
Rx = [1, 0, 0; 0, cos(roll), -sin(roll); 0, sin(roll), cos(roll)];
Ry = [cos(pitch), 0, sin(pitch); 0, 1, 0; -sin(pitch), 0, cos(pitch)];
Rz = [cos(yaw), -sin(yaw), 0; sin(yaw), cos(yaw), 0; 0, 0, 1];
R = Rz * Ry * Rx;
end

function a = steeringReference(directionCosines, Nx, Ny, spacingOverLambda)
vx = exp(1i * 2 * pi * spacingOverLambda(1) * directionCosines(1, :) .* (0:Nx-1).');
vy = exp(1i * 2 * pi * spacingOverLambda(2) * directionCosines(2, :) .* (0:Ny-1).');
a = zeros(Nx * Ny, size(directionCosines, 2));
for index = 1:size(directionCosines, 2)
    a(:, index) = kron(vx(:, index), vy(:, index));
end
end
