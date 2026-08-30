function validateGuLeshem2012(cfg)
%VALIDATEGULESHEM2012 校验Gu与Leshem（2012）复现及UPA迁移配置。
%   该校验器与主配置校验分文件保存，避免论文专用参数掺入通用场景
%   校验逻辑。未经修改的Rec需要直接求样本协方差逆，因此快拍数必须
%   不小于阵元数；这一限制用于阻止代码悄悄加入论文中没有的加载。

if ~isstruct(cfg.literature) || ...
        ~isfield(cfg.literature, 'guLeshem2012')
    experimentError('literature.guLeshem2012 is required.');
end
paper = cfg.literature.guLeshem2012;
required = {'sensorCount', 'spacingOverLambda', 'nominalDesiredDeg', ...
    'nominalInterfererDeg', 'interfererInrDb', 'desiredSectorDeg', ...
    'lookMismatchHalfWidthDeg', 'noisePower', 'fixedSnrDb', ...
    'fixedSnapshotCount', 'snrSweepDb', 'snapshotSweep', ...
    'monteCarloTrials', 'angleGridStepDeg', 'methodNames', ...
    'wcoSquaredRadiusScale', 'espSignalSubspaceDimension', ...
    'upaTransferSnapshotCount', 'upaTransferJitterDeg', ...
    'upaTransferGridSizes', 'upaTransferPeriods', 'upaDiagnosis', ...
    'experimentIndex'};
if ~isstruct(paper) || ~all(isfield(paper, required))
    experimentError('literature.guLeshem2012 is incomplete.');
end

validatePositiveIntegerCandidates(paper.sensorCount, ...
    'guLeshem2012.sensorCount');
if ~isscalar(paper.spacingOverLambda) || ...
        ~isfinite(paper.spacingOverLambda) || paper.spacingOverLambda <= 0
    experimentError('guLeshem2012.spacingOverLambda must be positive.');
end
angles = [paper.nominalDesiredDeg, paper.nominalInterfererDeg];
if ~isnumeric(angles) || numel(angles) ~= 3 || ...
        any(~isfinite(angles)) || any(abs(angles) > 90)
    experimentError('Gu-Leshem nominal angles must contain three valid ULA DOAs.');
end
if ~isnumeric(paper.nominalInterfererDeg) || ...
        ~isequal(size(paper.nominalInterfererDeg), [1, 2]) || ...
        ~isnumeric(paper.interfererInrDb) || ...
        ~isequal(size(paper.interfererInrDb), [1, 2]) || ...
        any(~isfinite(paper.interfererInrDb))
    experimentError('Gu-Leshem interferer directions and INR must be 1-by-2 vectors.');
end

sector = paper.desiredSectorDeg;
if ~isnumeric(sector) || ~isequal(size(sector), [1, 2]) || ...
        any(~isfinite(sector)) || sector(1) < -90 || sector(2) > 90 || ...
        sector(1) >= sector(2) || paper.nominalDesiredDeg < sector(1) || ...
        paper.nominalDesiredDeg > sector(2)
    experimentError('Gu-Leshem desired sector must contain the nominal desired DOA.');
end
if ~isscalar(paper.lookMismatchHalfWidthDeg) || ...
        ~isfinite(paper.lookMismatchHalfWidthDeg) || ...
        paper.lookMismatchHalfWidthDeg < 0
    experimentError('Gu-Leshem mismatch bound must be finite and nonnegative.');
end
if ~isscalar(paper.noisePower) || ~isfinite(paper.noisePower) || ...
        paper.noisePower <= 0 || ~isscalar(paper.fixedSnrDb) || ...
        ~isfinite(paper.fixedSnrDb)
    experimentError('Gu-Leshem noise power and fixed SNR are invalid.');
end

validatePositiveIntegerCandidates(paper.fixedSnapshotCount, ...
    'guLeshem2012.fixedSnapshotCount');
validatePositiveIntegerCandidates(paper.snapshotSweep, ...
    'guLeshem2012.snapshotSweep');
if paper.fixedSnapshotCount < paper.sensorCount || ...
        any(paper.snapshotSweep < paper.sensorCount)
    experimentError(['Gu-Leshem unmodified Rec requires every paper ' ...
        'snapshot count to be at least the ULA sensor count.']);
end
if ~isnumeric(paper.snrSweepDb) || ~isvector(paper.snrSweepDb) || ...
        isempty(paper.snrSweepDb) || any(~isfinite(paper.snrSweepDb))
    experimentError('guLeshem2012.snrSweepDb must be a finite nonempty vector.');
end
validatePositiveIntegerCandidates(paper.monteCarloTrials, ...
    'guLeshem2012.monteCarloTrials');
if ~isscalar(paper.angleGridStepDeg) || ...
        ~isfinite(paper.angleGridStepDeg) || paper.angleGridStepDeg <= 0 || ...
        abs((sector(1) + 90) / paper.angleGridStepDeg - ...
        round((sector(1) + 90) / paper.angleGridStepDeg)) > 1e-10 || ...
        abs((90 - sector(2)) / paper.angleGridStepDeg - ...
        round((90 - sector(2)) / paper.angleGridStepDeg)) > 1e-10
    experimentError('Gu-Leshem angle step must exactly divide both complement sectors.');
end

allowedMethods = ["oracle", "esp", "wco", "gu_leshem2012_rec"];
methods = reshape(string(paper.methodNames), 1, []);
if ~isstring(paper.methodNames) || isempty(methods) || ...
        numel(unique(methods)) ~= numel(methods) || ...
        any(~ismember(methods, allowedMethods))
    experimentError(['guLeshem2012.methodNames must be a unique nonempty ' ...
        'subset of Oracle, ESP, WCO, and Rec.']);
end
if ~ismember("gu_leshem2012_rec", methods)
    experimentError(['guLeshem2012.methodNames must include ' ...
        'gu_leshem2012_rec because both experiments diagnose this reconstruction.']);
end
if ~isscalar(paper.wcoSquaredRadiusScale) || ...
        ~isfinite(paper.wcoSquaredRadiusScale) || ...
        paper.wcoSquaredRadiusScale <= 0 || paper.wcoSquaredRadiusScale >= 1
    experimentError('Gu-Leshem WCO squared-radius scale must lie in (0,1).');
end
validatePositiveIntegerCandidates(paper.espSignalSubspaceDimension, ...
    'guLeshem2012.espSignalSubspaceDimension');
if paper.espSignalSubspaceDimension > paper.sensorCount
    experimentError('Gu-Leshem ESP dimension cannot exceed the ULA sensor count.');
end

validatePositiveIntegerCandidates(paper.upaTransferSnapshotCount, ...
    'guLeshem2012.upaTransferSnapshotCount');
if paper.upaTransferSnapshotCount < cfg.array.Nx * cfg.array.Ny
    experimentError('UPA transfer requires at least M snapshots for unmodified Rec.');
end
validateNonnegativeCandidates(paper.upaTransferJitterDeg, ...
    'guLeshem2012.upaTransferJitterDeg', inf);
validatePositiveIntegerCandidates(paper.upaTransferGridSizes, ...
    'guLeshem2012.upaTransferGridSizes');
if any(paper.upaTransferGridSizes < 3) || ...
        any(mod(paper.upaTransferGridSizes, 2) == 0)
    experimentError('UPA transfer grid sizes must be odd integers at least three.');
end
validatePositiveIntegerCandidates(paper.upaTransferPeriods, ...
    'guLeshem2012.upaTransferPeriods');
if paper.upaTransferPeriods > cfg.timeline.K
    experimentError('UPA transfer periods cannot exceed timeline.K.');
end
validatePositiveIntegerCandidates(paper.experimentIndex, ...
    'guLeshem2012.experimentIndex');
validateUpaDiagnosis(paper.upaDiagnosis, cfg);
end

function validateUpaDiagnosis(diagnosis, cfg)
%VALIDATEUPADIAGNOSIS 校验双快拍UPA积分诊断的固定参数合同。
required = {'snapshotCounts', 'gridSizes', 'measureNames', ...
    'selectedGridSize', 'selectedMeasureName', 'loadingRatio', ...
    'fitRegularization', 'peakThresholdRatio', 'peakNmsRadiusBins', ...
    'maximumCandidateCount', 'jitterSweepDeg', 'monteCarloTrials', ...
    'periods', 'experimentIndex'};
if ~isstruct(diagnosis) || ~all(isfield(diagnosis, required))
    experimentError('guLeshem2012.upaDiagnosis is incomplete.');
end

if ~isequal(diagnosis.snapshotCounts, [8, 64])
    experimentError('UPA diagnosis snapshot counts must equal [8, 64].');
end
validatePositiveIntegerCandidates(diagnosis.gridSizes, ...
    'guLeshem2012.upaDiagnosis.gridSizes');
if any(diagnosis.gridSizes < 3) || any(mod(diagnosis.gridSizes, 2) == 0)
    experimentError('UPA diagnosis grid sizes must be odd integers at least three.');
end
validatePositiveIntegerScalar(diagnosis.selectedGridSize, ...
    'guLeshem2012.upaDiagnosis.selectedGridSize');
if ~ismember(diagnosis.selectedGridSize, diagnosis.gridSizes)
    experimentError('UPA diagnosis selected grid size must be a configured grid size.');
end

allowedMeasures = ["direction_cosine", "solid_angle", "az_el_solid_angle"];
measures = reshape(string(diagnosis.measureNames), 1, []);
if ~isstring(diagnosis.measureNames) || isempty(measures) || ...
        numel(unique(measures)) ~= numel(measures) || ...
        any(~ismember(measures, allowedMeasures))
    experimentError(['UPA diagnosis measure names must be a nonempty unique ' ...
        'subset of the three supported measures.']);
end
if ~isstring(diagnosis.selectedMeasureName) || ...
        ~isscalar(diagnosis.selectedMeasureName) || ...
        ~ismember(diagnosis.selectedMeasureName, measures)
    experimentError('UPA diagnosis selected measure must be a configured measure name.');
end

validatePositiveFiniteScalar(diagnosis.loadingRatio, ...
    'guLeshem2012.upaDiagnosis.loadingRatio');
validateNonnegativeFiniteScalar(diagnosis.fitRegularization, ...
    'guLeshem2012.upaDiagnosis.fitRegularization');
validateUnitIntervalScalar(diagnosis.peakThresholdRatio, ...
    'guLeshem2012.upaDiagnosis.peakThresholdRatio');
validateNonnegativeIntegerScalar(diagnosis.peakNmsRadiusBins, ...
    'guLeshem2012.upaDiagnosis.peakNmsRadiusBins');
validatePositiveIntegerScalar(diagnosis.maximumCandidateCount, ...
    'guLeshem2012.upaDiagnosis.maximumCandidateCount');
validateNonnegativeCandidates(diagnosis.jitterSweepDeg, ...
    'guLeshem2012.upaDiagnosis.jitterSweepDeg', inf);
validatePositiveIntegerScalar(diagnosis.monteCarloTrials, ...
    'guLeshem2012.upaDiagnosis.monteCarloTrials');
validatePositiveIntegerScalar(diagnosis.periods, ...
    'guLeshem2012.upaDiagnosis.periods');
if diagnosis.periods > cfg.timeline.K
    experimentError('UPA diagnosis periods cannot exceed timeline.K.');
end
validatePositiveIntegerScalar(diagnosis.experimentIndex, ...
    'guLeshem2012.upaDiagnosis.experimentIndex');
end

function validatePositiveIntegerCandidates(values, name)
%VALIDATEPOSITIVEINTEGERCANDIDATES 校验非空、唯一的正整数向量。
if ~isnumeric(values) || ~isvector(values) || isempty(values) || ...
        any(~isfinite(values)) || any(values <= 0) || ...
        any(values ~= floor(values)) || numel(unique(values)) ~= numel(values)
    experimentError('%s must contain unique positive integers.', name);
end
end

function validateNonnegativeCandidates(values, name, exclusiveUpperBound)
%VALIDATENONNEGATIVECANDIDATES 校验非空、唯一且落在规定范围内的向量。
if ~isnumeric(values) || ~isvector(values) || isempty(values) || ...
        any(~isfinite(values)) || any(values < 0) || ...
        any(values >= exclusiveUpperBound) || ...
        numel(unique(values)) ~= numel(values)
    experimentError('%s contains invalid candidates.', name);
end
end

function validatePositiveFiniteScalar(value, name)
%VALIDATEPOSITIVEFINITESCALAR 校验有限正标量。
if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value) || value <= 0
    experimentError('%s must be a positive finite scalar.', name);
end
end

function validateNonnegativeFiniteScalar(value, name)
%VALIDATENONNEGATIVEFINITESCALAR 校验有限非负标量。
if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value) || value < 0
    experimentError('%s must be a nonnegative finite scalar.', name);
end
end

function validateUnitIntervalScalar(value, name)
%VALIDATEUNITINTERVALSCALAR 校验闭区间[0,1]内的有限标量。
if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value) || ...
        value < 0 || value > 1
    experimentError('%s must be a finite scalar in [0, 1].', name);
end
end

function validatePositiveIntegerScalar(value, name)
%VALIDATEPOSITIVEINTEGERSCALAR 校验有限正整数标量。
if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value) || ...
        value <= 0 || value ~= floor(value)
    experimentError('%s must be a positive integer scalar.', name);
end
end

function validateNonnegativeIntegerScalar(value, name)
%VALIDATENONNEGATIVEINTEGERSCALAR 校验有限非负整数标量。
if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value) || ...
        value < 0 || value ~= floor(value)
    experimentError('%s must be a nonnegative integer scalar.', name);
end
end

function experimentError(messageText, varargin)
%EXPERIMENTERROR 使用论文实验配置专用的统一错误标识。
error('uav:config:InvalidExperimentParameter', messageText, varargin{:});
end
