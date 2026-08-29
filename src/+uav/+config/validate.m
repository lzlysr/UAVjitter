function cfg = validate(cfg)
%VALIDATE 校验配置并计算阵列、场景、功率和频谱网格派生量。
%   本函数既是输入防线，也是“原始参数->模型可直接使用参数”的转换层。
%   它不修改用户给出的物理含义，只规范空干扰向量和输出路径，并生成
%   cfg.derived。重复调用必须得到一致结果。

requiredSections = {'random', 'array', 'timeline', 'attitude', 'sources', ...
    'signal', 'spectrum', 'incm', 'beamforming', 'monteCarlo', 'experiments', ...
    'plot', 'output'};
for index = 1:numel(requiredSections)
    if ~isfield(cfg, requiredSections{index})
        error('uav:config:MissingSection', 'Missing configuration section "%s".', ...
            requiredSections{index});
    end
end

% 第一组：阵列、时间、姿态、随机和网格的基本尺寸与范围检查。
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
% 第二组：功率、实验扫描和结果管理参数检查。
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

% 第三组：方向采用[azimuth,elevation]，当前下行场景限定在阵列下半球。
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
% J由干扰方向行数自动推导；J=0时统一保存为1×0功率向量。
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
% 名义姿态下预计算世界方向、阵列方向、方向余弦和所有源的名义流形。
worldUnitVectors = azElToUnitVectors(azElDeg);
nominalRotation = rotationZYX(deg2rad(cfg.attitude.nominalDeg));
nominalArrayUnitVectors = nominalRotation.' * worldUnitVectors;
directionCosines = nominalArrayUnitVectors(1:2, :);
nominalSteering = steeringReference(directionCosines, cfg.array.Nx, cfg.array.Ny, ...
    cfg.array.spacingOverLambda);
[gridAxis, xiX, xiY, visibleMask, visibleLinearIndices, ...
    gridDirectionCosines, gridManifold, gridCellArea] = spectrumGrid(cfg);

% SNR/INR均相对于配置噪声功率定义，进入模型后转换为线性功率。
noisePower = cfg.signal.noisePower;
desiredPower = noisePower * 10^(cfg.signal.snrDb / 10);
interfererPowers = noisePower * 10.^(cfg.sources.interfererInrDb(:).' / 10);
defaultRootDir = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
% 若调用者显式配置临时输出目录，则保留该目录；否则落到工程根目录。
rootDir = configuredOutputPath(cfg.output, 'rootDir', defaultRootDir);
cfg.output = struct( ...
    'rootDir', rootDir, ...
    'cacheDir', configuredOutputPath(cfg.output, 'cacheDir', ...
        fullfile(rootDir, 'data', 'cache')), ...
    'rawDir', configuredOutputPath(cfg.output, 'rawDir', ...
        fullfile(rootDir, 'results', 'raw')), ...
    'figureDir', configuredOutputPath(cfg.output, 'figureDir', ...
        fullfile(rootDir, 'results', 'figures')), ...
    'tableDir', configuredOutputPath(cfg.output, 'tableDir', ...
        fullfile(rootDir, 'results', 'tables')), ...
    'logDir', configuredOutputPath(cfg.output, 'logDir', ...
        fullfile(rootDir, 'results', 'logs')), ...
    'overwrite', cfg.output.overwrite);
% 派生量只由上面的已校验原始参数计算，不允许算法反向修改。
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
    'nominalSteeringReference', nominalSteering(:, 1), ...
    'spectrumGrid', struct( ...
        'axis', gridAxis, ...
        'xiX', xiX, ...
        'xiY', xiY, ...
        'visibleMask', visibleMask, ...
        'visibleLinearIndices', visibleLinearIndices, ...
        'directionCosines', gridDirectionCosines, ...
        'manifold', gridManifold, ...
        'cellArea', gridCellArea));
end

function [axisValues, xiX, xiY, visibleMask, visibleIndices, ...
        directionCosines, manifold, cellArea] = spectrumGrid(cfg)
%SPECTRUMGRID 预计算二维方向余弦网格和所有可见节点的UPA流形。
% 缓存M×Nvisible流形可避免每个周期、每种方法重复生成相同导向矢量。
axisValues = linspace(-1, 1, cfg.spectrum.gridSize);
[xiX, xiY] = meshgrid(axisValues, axisValues);
visibleMask = xiX.^2 + xiY.^2 <= 1;
visibleIndices = find(visibleMask);
directionCosines = [xiX(visibleIndices).'; xiY(visibleIndices).'];
manifold = steeringReference(directionCosines, cfg.array.Nx, cfg.array.Ny, ...
    cfg.array.spacingOverLambda);
gridStep = axisValues(2) - axisValues(1);
cellArea = gridStep^2;
end

function pathValue = configuredOutputPath(output, field, defaultValue)
%CONFIGUREDOUTPUTPATH 保留合法的调用者自定义目录，否则使用工程默认目录。
if ~isfield(output, field) || isempty(output.(field))
    pathValue = defaultValue;
    return;
end
candidate = output.(field);
if ~(ischar(candidate) && isrow(candidate)) && ...
        ~(isstring(candidate) && isscalar(candidate) && ~ismissing(candidate))
    error('uav:config:InvalidExperimentParameter', ...
        'output.%s must be a nonmissing character vector or string scalar.', field);
end
pathValue = char(candidate);
if isempty(strtrim(pathValue))
    error('uav:config:InvalidExperimentParameter', ...
        'output.%s must not be empty.', field);
end
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
        ~islogical(cfg.plot.closeAfterExport) || ...
        ~isscalar(cfg.plot.closeAfterExport) || ...
        ~islogical(cfg.plot.saveFigures) || ~isscalar(cfg.plot.saveFigures) || ...
        ~cfg.plot.saveFigures || ...
        ~isstring(cfg.plot.formats) || ...
        ~isequal(reshape(cfg.plot.formats, 1, []), ["fig", "pdf", "png"])
    experimentError('Plot switches and formats are invalid.');
end
positivePlotFields = {'resolutionDpi', 'lineWidth', 'markerSize', ...
    'fontSize', 'maxTraceTrials'};
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
