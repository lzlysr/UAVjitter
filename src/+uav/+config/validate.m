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
validateExperimentParameters(cfg);

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
    'logDir', fullfile(rootDir, 'results', 'logs'), ...
    'overwrite', cfg.output.overwrite);
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
validateAlgorithmIntegerVector(cfg.spectrum.alignmentMaxShiftBins, 2, ...
    'spectrum.alignmentMaxShiftBins', false);
validateUnitIntervalScalar(cfg.spectrum.alignmentMinCorrelation, ...
    'spectrum.alignmentMinCorrelation', true);
validateUnitIntervalScalar(cfg.spectrum.alignmentAmbiguityMargin, ...
    'spectrum.alignmentAmbiguityMargin', true);
validateAlgorithmIntegerScalar(cfg.spectrum.alignmentMinimumOverlap, ...
    'spectrum.alignmentMinimumOverlap', true);
if cfg.spectrum.alignmentMinimumOverlap < 2
    error('uav:config:InvalidAlgorithmParameter', ...
        'spectrum.alignmentMinimumOverlap must be at least two.');
end
validateUnitIntervalScalar(cfg.spectrum.peakThresholdRatio, ...
    'spectrum.peakThresholdRatio', true);
validateAlgorithmIntegerScalar(cfg.spectrum.peakNmsRadiusBins, ...
    'spectrum.peakNmsRadiusBins', false);
validateAlgorithmIntegerScalar(cfg.spectrum.quadratureCoarseStride, ...
    'spectrum.quadratureCoarseStride', true);
validateAlgorithmIntegerScalar(cfg.spectrum.quadratureFineRadiusBins, ...
    'spectrum.quadratureFineRadiusBins', false);
validateAlgorithmIntegerScalar(cfg.spectrum.quadratureRefinement, ...
    'spectrum.quadratureRefinement', true);
if cfg.spectrum.quadratureRefinement <= 1
    error('uav:config:InvalidAlgorithmParameter', ...
        'spectrum.quadratureRefinement must be greater than one.');
end
validatePositiveFiniteScalar(cfg.spectrum.quadratureDedupTolerance, ...
    'spectrum.quadratureDedupTolerance');
gridStep = 2 / (cfg.spectrum.gridSize - 1);
maximumDedupTolerance = gridStep / (2 * cfg.spectrum.quadratureRefinement);
if cfg.spectrum.quadratureDedupTolerance > maximumDedupTolerance
    error('uav:config:InvalidAlgorithmParameter', ...
        ['spectrum.quadratureDedupTolerance must not exceed half the ' ...
        'fine-node spacing.']);
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
validateAblationSwitches(cfg.incm.ablation);
validateNonnegativeFiniteScalar(cfg.beamforming.diagonalLoading, ...
    'beamforming.diagonalLoading');
validateNonnegativeFiniteScalar(cfg.beamforming.epsilonS, 'beamforming.epsilonS');
validateNonnegativeFiniteScalar(cfg.beamforming.gamma, 'beamforming.gamma');
validatePositiveFiniteScalar(cfg.beamforming.covariancePsdTolerance, ...
    'beamforming.covariancePsdTolerance');
validatePositiveFiniteScalar(cfg.beamforming.minimumReciprocalCondition, ...
    'beamforming.minimumReciprocalCondition');
if cfg.beamforming.minimumReciprocalCondition >= 1
    error('uav:config:InvalidAlgorithmParameter', ...
        'beamforming.minimumReciprocalCondition must lie in (0, 1).');
end
validatePositiveFiniteScalar(cfg.beamforming.solver.optimalityTolerance, ...
    'beamforming.solver.optimalityTolerance');
validatePositiveFiniteScalar(cfg.beamforming.solver.constraintTolerance, ...
    'beamforming.solver.constraintTolerance');
validatePositiveFiniteScalar(cfg.beamforming.solver.stepTolerance, ...
    'beamforming.solver.stepTolerance');
validatePositiveFiniteScalar(cfg.beamforming.solver.initialInteriorMargin, ...
    'beamforming.solver.initialInteriorMargin');
validatePositiveFiniteScalar(cfg.beamforming.solver.feasibilityTolerance, ...
    'beamforming.solver.feasibilityTolerance');
validateAlgorithmIntegerScalar(cfg.beamforming.solver.maxIterations, ...
    'beamforming.solver.maxIterations', true);
validateFiniteVector(cfg.beamforming.calibration.designJitterStdDeg, 3, ...
    'beamforming.calibration.designJitterStdDeg');
if any(cfg.beamforming.calibration.designJitterStdDeg < 0)
    error('uav:config:InvalidAlgorithmParameter', ...
        'beamforming.calibration.designJitterStdDeg must be nonnegative.');
end
validateAlgorithmFiniteScalar(cfg.beamforming.calibration.quantile, ...
    'beamforming.calibration.quantile');
if cfg.beamforming.calibration.quantile <= 0 || ...
        cfg.beamforming.calibration.quantile >= 1
    error('uav:config:InvalidAlgorithmParameter', ...
        'beamforming.calibration.quantile must lie in (0, 1).');
end
validateAlgorithmIntegerScalar(cfg.beamforming.calibration.sampleCount, ...
    'beamforming.calibration.sampleCount', true);
validateAlgorithmIntegerScalar(cfg.beamforming.calibration.substream, ...
    'beamforming.calibration.substream', true);
end

function validateExperimentParameters(cfg)
expectedMethods = ["oracle", "frozen", "scm", "dl", "static", "temporal"];
expectedAblations = ["full", "lambda_t_zero", "identity_alignment", ...
    "uniform_quadrature", "no_structure", "pilot_projection"];
if ~isstring(cfg.experiments.methodNames) || ...
        ~isequal(reshape(cfg.experiments.methodNames, 1, []), expectedMethods) || ...
        ~isstring(cfg.experiments.ablationNames) || ...
        ~isequal(reshape(cfg.experiments.ablationNames, 1, []), expectedAblations)
    experimentError('Method and ablation names must match the fixed Task 5 contract.');
end

periodFields = {'modelValidation', 'jitterSweep', 'snapshotSweep', ...
    'dynamicStability', 'ablation'};
if ~isstruct(cfg.experiments.periodCounts) || ...
        ~all(isfield(cfg.experiments.periodCounts, periodFields))
    experimentError('experiments.periodCounts is incomplete.');
end
for index = 1:numel(periodFields)
    value = cfg.experiments.periodCounts.(periodFields{index});
    if ~isscalar(value) || ~isnumeric(value) || ~isfinite(value) || ...
            value <= 0 || value ~= floor(value) || value > cfg.timeline.K
        experimentError('Every experiment period count must be a positive integer no larger than K.');
    end
end
if ~isscalar(cfg.experiments.bootstrapReplicates) || ...
        ~isfinite(cfg.experiments.bootstrapReplicates) || ...
        cfg.experiments.bootstrapReplicates <= 0 || ...
        cfg.experiments.bootstrapReplicates ~= floor(cfg.experiments.bootstrapReplicates)
    experimentError('Experiment replicate counts must be positive integers.');
end
if ~isscalar(cfg.experiments.outageThresholdDb) || ...
        ~isfinite(cfg.experiments.outageThresholdDb)
    experimentError('The outage threshold must be a finite dB scalar.');
end
if ~isscalar(cfg.experiments.oracleTolerance) || ...
        ~isfinite(cfg.experiments.oracleTolerance) || ...
        cfg.experiments.oracleTolerance < 0
    experimentError('The Oracle comparison tolerance must be nonnegative and finite.');
end
aggregation = cfg.experiments.aggregation;
if ~isstruct(aggregation) || ~all(isfield(aggregation, ...
        {'medianQuantile', 'lowQuantile', 'confidenceQuantiles', 'runtimeQuantile'})) || ...
        aggregation.medianQuantile ~= 0.5 || ...
        ~isscalar(aggregation.lowQuantile) || aggregation.lowQuantile <= 0 || ...
        aggregation.lowQuantile >= aggregation.medianQuantile || ...
        ~isnumeric(aggregation.confidenceQuantiles) || ...
        ~isequal(size(aggregation.confidenceQuantiles), [1, 2]) || ...
        any(~isfinite(aggregation.confidenceQuantiles)) || ...
        aggregation.confidenceQuantiles(1) <= 0 || ...
        aggregation.confidenceQuantiles(2) >= 1 || ...
        aggregation.confidenceQuantiles(1) >= aggregation.confidenceQuantiles(2) || ...
        ~isscalar(aggregation.runtimeQuantile) || ...
        aggregation.runtimeQuantile <= 0 || aggregation.runtimeQuantile >= 1
    experimentError('Experiment aggregation controls are invalid.');
end

substreams = cfg.experiments.substreams;
substreamFields = {'experimentStride', 'pointStride', 'trialStride', ...
    'scenarioOffset', 'snapshotOffset', 'bootstrapOffset'};
if ~isstruct(substreams) || ~all(isfield(substreams, substreamFields))
    experimentError('The substream allocation scheme is incomplete.');
end
for index = 1:numel(substreamFields)
    value = substreams.(substreamFields{index});
    if ~isscalar(value) || ~isnumeric(value) || ~isfinite(value) || ...
            value <= 0 || value ~= floor(value)
        experimentError('Substream controls must be positive integers.');
    end
end
maximumExperimentPeriods = max(structfun(@double, cfg.experiments.periodCounts));
if substreams.snapshotOffset <= substreams.scenarioOffset || ...
        substreams.bootstrapOffset <= substreams.snapshotOffset + maximumExperimentPeriods || ...
        substreams.trialStride <= substreams.bootstrapOffset || ...
        substreams.pointStride <= cfg.monteCarlo.nTrials * substreams.trialStride || ...
        substreams.experimentStride <= ...
            max(numel(cfg.experiments.jitterSweepDeg), ...
            numel(cfg.experiments.snapshotSweep)) * substreams.pointStride
    experimentError('Substream strides do not isolate experiments, points, trials, and periods.');
end
if ~isstring(cfg.experiments.codeVersion) || ~isscalar(cfg.experiments.codeVersion) || ...
        strlength(cfg.experiments.codeVersion) == 0
    experimentError('The host-independent code version identifier is required.');
end
if ~islogical(cfg.plot.visible) || ~isscalar(cfg.plot.visible) || ...
        ~islogical(cfg.plot.saveFigures) || ~isscalar(cfg.plot.saveFigures) || ...
        ~isstring(cfg.plot.formats) || ...
        ~isequal(reshape(cfg.plot.formats, 1, []), ["fig", "pdf", "png"])
    experimentError('Plot switches and formats are invalid.');
end
positivePlotFields = {'resolutionDpi', 'lineWidth', 'markerSize', 'maxTraceTrials'};
for index = 1:numel(positivePlotFields)
    value = cfg.plot.(positivePlotFields{index});
    if ~isscalar(value) || ~isnumeric(value) || ~isfinite(value) || value <= 0
        experimentError('Plot numeric controls must be positive finite scalars.');
    end
end
if cfg.plot.maxTraceTrials ~= floor(cfg.plot.maxTraceTrials) || ...
        ~isnumeric(cfg.plot.figurePosition) || numel(cfg.plot.figurePosition) ~= 4 || ...
        any(~isfinite(cfg.plot.figurePosition)) || any(cfg.plot.figurePosition(3:4) <= 0)
    experimentError('Plot trace count or figure position is invalid.');
end
if ~isfield(cfg.output, 'overwrite') || ~islogical(cfg.output.overwrite) || ...
        ~isscalar(cfg.output.overwrite)
    experimentError('output.overwrite must be a logical scalar.');
end
end

function experimentError(messageText)
error('uav:config:InvalidExperimentParameter', '%s', messageText);
end

function validateAblationSwitches(ablation)
required = {'temporalContinuity', 'alignment', 'adaptiveQuadrature', ...
    'structureEnhancement', 'pilotProjection'};
if ~isstruct(ablation) || ~all(isfield(ablation, required))
    error('uav:config:InvalidAlgorithmParameter', ...
        'incm.ablation must provide every Task 4 switch.');
end
for index = 1:numel(required)
    value = ablation.(required{index});
    if ~islogical(value) || ~isscalar(value)
        error('uav:config:InvalidAlgorithmParameter', ...
            'incm.ablation.%s must be a logical scalar.', required{index});
    end
end
end

function validateUnitIntervalScalar(value, name, allowZero)
validateAlgorithmFiniteScalar(value, name);
if value > 1 || (allowZero && value < 0) || (~allowZero && value <= 0)
    error('uav:config:InvalidAlgorithmParameter', ...
        '%s must lie in the configured unit interval.', name);
end
end

function validateAlgorithmIntegerVector(value, expectedLength, name, positive)
if ~isnumeric(value) || ~isvector(value) || numel(value) ~= expectedLength || ...
        any(~isfinite(value)) || any(value ~= floor(value)) || ...
        (positive && any(value <= 0)) || (~positive && any(value < 0))
    error('uav:config:InvalidAlgorithmParameter', ...
        '%s must be an integer vector with the configured sign.', name);
end
end

function validateAlgorithmIntegerScalar(value, name, positive)
if ~isscalar(value) || ~isnumeric(value) || ~isfinite(value) || ...
        value ~= floor(value) || (positive && value <= 0) || (~positive && value < 0)
    error('uav:config:InvalidAlgorithmParameter', ...
        '%s must be an integer with the configured sign.', name);
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
