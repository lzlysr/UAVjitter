classdef TestTask5Experiments < matlab.unittest.TestCase
    methods (TestMethodSetup)
        function addRuntimePaths(testCase)
            rootDir = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(rootDir));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(rootDir, 'src')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(rootDir, 'experiments')));
        end
    end

    methods (Test)
        function configDefinesAndValidatesAllExperimentControls(testCase)
            cfg = config('smoke');

            testCase.verifyEqual(cfg.experiments.methodNames, ...
                ["oracle", "frozen", "scm", "dl", "static", "temporal"]);
            testCase.verifyEqual(cfg.experiments.ablationNames, ...
                ["full", "lambda_t_zero", "identity_alignment", ...
                "uniform_quadrature", "no_structure", "pilot_projection"]);
            testCase.verifyEqual(cfg.experiments.periodCounts.dynamicStability, ...
                cfg.timeline.K);
            testCase.verifyGreaterThan(cfg.experiments.bootstrapReplicates, 0);
            testCase.verifyTrue(isfinite(cfg.experiments.outageThresholdDb));
            testCase.verifyGreaterThan(cfg.experiments.substreams.experimentStride, ...
                cfg.experiments.substreams.pointStride);
            testCase.verifyEqual(cfg.plot.formats, ["fig", "pdf", "png"]);
            testCase.verifyTrue(cfg.output.overwrite);

            invalid = cfg;
            invalid.experiments.methodNames(2) = "oracle";
            testCase.verifyError(@() uav.config.validate(invalid), ...
                'uav:config:InvalidExperimentParameter');
            invalid = cfg;
            invalid.experiments.periodCounts.jitterSweep = 0;
            testCase.verifyError(@() uav.config.validate(invalid), ...
                'uav:config:InvalidExperimentParameter');
            invalid = cfg;
            invalid.experiments.substreams.trialStride = ...
                invalid.experiments.substreams.snapshotOffset;
            testCase.verifyError(@() uav.config.validate(invalid), ...
                'uav:config:InvalidExperimentParameter');
            invalid = cfg;
            invalid.plot.saveFigures = false;
            testCase.verifyError(@() uav.config.validate(invalid), ...
                'uav:config:InvalidExperimentParameter');
        end

        function linearAggregationAndBootstrapAreDeterministic(testCase)
            cfg = compactConfig(testCase);
            values = [1, 2; 4, NaN; 9, 8; 16, 32];

            first = uav.experiment.summarizeLinear(values, cfg, 71);
            second = uav.experiment.summarizeLinear(values, cfg, 71);

            testCase.verifyEqual(first, second);
            testCase.verifyEqual(first.medianLinear, [6.5, 8], 'AbsTol', 0);
            testCase.verifyEqual(first.medianDb(1), 10 * log10(6.5), 'AbsTol', 1e-12);
            testCase.verifyNotEqual(first.medianDb(1), median(10 * log10(values(:, 1))));
            testCase.verifySize(first.ci95Linear, [2, 2]);
            testCase.verifySize(first.ci95Db, [2, 2]);
            testCase.verifyEqual(first.validCounts, [4, 3]);
        end

        function trialClusterBootstrapPreservesMethodPairing(testCase)
            cfg = compactConfig(testCase);
            base = reshape([1, 2, 3, 4, 8, 16], 3, 2);
            values = cat(3, base, 4 * base);

            summary = uav.experiment.summarizeTrials(values, cfg, 73);

            testCase.verifyEqual(summary.bootstrapUnit, "trial");
            testCase.verifyEqual(summary.independentTrialCount, 3);
            testCase.verifyEqual(summary.periodsPerTrial, 2);
            testCase.verifyEqual(summary.ci95Linear(:, 2), ...
                4 * summary.ci95Linear(:, 1), 'AbsTol', 0);
            testCase.verifyEqual(summary.ci95Db(:, 2) - ...
                summary.ci95Db(:, 1), 10 * log10(4) * ones(2, 1), ...
                'AbsTol', 1e-12);
        end

        function paperSizedAggregationDoesNotBootstrapIndividualPeriods(testCase)
            cfg = config('paper');
            cfg.experiments.bootstrapReplicates = 2;
            cfg = uav.config.validate(cfg);
            values = ones(cfg.monteCarlo.nTrials, cfg.timeline.K, 2);
            values(:, :, 2) = 2;

            summary = uav.experiment.summarizeTrials(values, cfg, 75);

            testCase.verifyEqual(summary.independentTrialCount, 500);
            testCase.verifyEqual(summary.periodsPerTrial, 500);
            testCase.verifyEqual(summary.ci95Linear, [1, 2; 1, 2], 'AbsTol', 0);
        end

        function sixMethodRunnerUsesOneSnapshotAndFixedCalibration(testCase)
            cfg = compactConfig(testCase);
            cfg.beamforming.epsilonS = 0.1;
            cfg = uav.config.validate(cfg);
            scenario = uav.scenario.generate(cfg, uav.util.makeStream(701, 1));
            obs1 = uav.scenario.generatePeriod(cfg, scenario, 1, ...
                uav.util.makeStream(701, 2));
            obs2 = uav.scenario.generatePeriod(cfg, scenario, 2, ...
                uav.util.makeStream(701, 3));

            [record1, state] = uav.experiment.runMethods(obs1, [], cfg);
            frozenWeight = state.frozenWeight;
            [record2, state] = uav.experiment.runMethods(obs2, state, cfg);

            testCase.verifyEqual(record1.methodNames, cfg.experiments.methodNames);
            testCase.verifySize(record1.sinrLinear, [1, 6]);
            testCase.verifySize(record1.runtimeSec, [1, 6]);
            testCase.verifyEqual(record1.snapshotSignatures, ...
                repmat(record1.snapshotSignatures(1, :), 6, 1), 'AbsTol', 0);
            testCase.verifyEqual(record2.snapshotSignatures, ...
                repmat(record2.snapshotSignatures(1, :), 6, 1), 'AbsTol', 0);
            testCase.verifyEqual(state.frozenWeight, frozenWeight, 'AbsTol', 0);
            testCase.verifyEqual(record1.epsilonS, 0.1, 'AbsTol', 0);
            testCase.verifyGreaterThanOrEqual(record1.sinrLinear(1), ...
                max(record1.sinrLinear(2:end)) - 1e-9);
            testCase.verifyTrue(all(isfinite(record1.sinrLinear)));
            testCase.verifyTrue(all(isfinite(record1.desiredGain)));
            testCase.verifyTrue(all(isfinite(record1.interferenceLeakage)));
        end

        function storedSubstreamsRegenerateEveryPairedSnapshot(testCase)
            cfg = compactConfig(testCase);
            cfg.monteCarlo.nTrials = 2;
            cfg.beamforming.epsilonS = 0.1;
            cfg = uav.config.validate(cfg);

            [raw, seeds] = uav.experiment.collectMethods(cfg, 2, 1, 2);

            for trialIndex = 1:cfg.monteCarlo.nTrials
                scenario = uav.scenario.generate(cfg, uav.util.makeStream( ...
                    cfg.random.baseSeed, seeds.scenarioSubstreams(trialIndex)));
                for periodIndex = 1:2
                    obs = uav.scenario.generatePeriod(cfg, scenario, periodIndex, ...
                        uav.util.makeStream(cfg.random.baseSeed, ...
                        seeds.snapshotSubstreams(trialIndex, periodIndex)));
                    expected = uav.util.snapshotSignature(obs.X);
                    stored = squeeze(raw.snapshotSignatures( ...
                        trialIndex, periodIndex, :, :));
                    testCase.verifyEqual(stored, repmat(expected, ...
                        numel(cfg.experiments.methodNames), 1), 'AbsTol', 0);
                end
            end
        end

        function algorithmObservationWhitelistsNonOracleInputs(testCase)
            cfg = compactConfig(testCase);
            scenario = uav.scenario.generate(cfg, uav.util.makeStream(719, 1));
            obs = uav.scenario.generatePeriod(cfg, scenario, 1, ...
                uav.util.makeStream(719, 2));
            obs.unapprovedDiagnostic = obs.truth.trueIncm;

            safe = uav.experiment.algorithmObservation(obs);

            testCase.verifyFalse(isfield(safe, 'truth'));
            testCase.verifyFalse(isfield(safe, 'unapprovedDiagnostic'));
            testCase.verifyEqual(sort(string(fieldnames(safe))), sort([ ...
                "X"; "Rx"; "pilot"; "k"; "nominalSteeringReference"]));
            testCase.verifyEqual(safe.X, obs.X, 'AbsTol', 0);
            testCase.verifyEqual(safe.Rx, obs.Rx, 'AbsTol', 0);
            testCase.verifyEqual(safe.pilot, obs.pilot, 'AbsTol', 0);
            testCase.verifyEqual(safe.nominalSteeringReference, ...
                obs.nominalSteeringReference, 'AbsTol', 0);
        end

        function ablationRunnerChangesExactlyOneConfiguredSwitch(testCase)
            cfg = compactConfig(testCase);
            cfg.beamforming.epsilonS = 0.1;
            cfg = uav.config.validate(cfg);
            scenario = uav.scenario.generate(cfg, uav.util.makeStream(709, 1));
            obs = uav.scenario.generatePeriod(cfg, scenario, 1, ...
                uav.util.makeStream(709, 2));

            [record, ~] = uav.experiment.runVariants(obs, [], cfg);

            testCase.verifyEqual(record.variantNames, cfg.experiments.ablationNames);
            testCase.verifyEqual(record.switchDeltaFields, ...
                ["none", "incm.lambdaT", "incm.ablation.alignment", ...
                "incm.ablation.adaptiveQuadrature", ...
                "incm.ablation.structureEnhancement", ...
                "incm.ablation.pilotProjection"]);
            testCase.verifyEqual(record.snapshotSignatures, ...
                repmat(record.snapshotSignatures(1, :), 6, 1), 'AbsTol', 0);
            testCase.verifyEqual(record.epsilonS, 0.1, 'AbsTol', 0);
            testCase.verifyTrue(all(isfinite(record.sinrLinear)));

            base = record.effectiveIncmConfigs{1};
            expected = base; expected.lambdaT = 0;
            testCase.verifyEqual(record.effectiveIncmConfigs{2}, expected);
            expected = base; expected.ablation.alignment = false;
            testCase.verifyEqual(record.effectiveIncmConfigs{3}, expected);
            expected = base; expected.ablation.adaptiveQuadrature = false;
            testCase.verifyEqual(record.effectiveIncmConfigs{4}, expected);
            expected = base; expected.ablation.structureEnhancement = false;
            testCase.verifyEqual(record.effectiveIncmConfigs{5}, expected);
            expected = base; expected.ablation.pilotProjection = true;
            testCase.verifyEqual(record.effectiveIncmConfigs{6}, expected);

            invalidCfg = cfg;
            invalidCfg.incm.lambdaT = 0;
            testCase.verifyError(@() uav.experiment.runVariants(obs, [], invalidCfg), ...
                'uav:experiment:InvalidFullAblationBaseline');
        end

        function metadataAndArtifactWriterCreateCompleteRerunnableBundle(testCase)
            cfg = compactConfig(testCase);
            fig = figure('Visible', 'off');
            plot([0, 1], [1, 2]);
            testCase.addTeardown(@() closeIfValid(fig));
            metadata = uav.experiment.metadata(cfg, "writer_contract", ...
                struct('trials', 1, 'periods', 2), 0.25);
            result = struct('cfg', cfg, 'metadata', metadata, ...
                'experimentName', "writer_contract", ...
                'methodNames', cfg.experiments.methodNames, ...
                'raw', struct('sinrLinear', ones(1, 2, 6)), ...
                'summary', struct('medianLinear', ones(1, 6)), ...
                'timestamp', metadata.generatedAt, 'runtimeSec', 0.25);
            summary = table((1:2).', [1; 2], 'VariableNames', {'index', 'value'});

            first = uav.experiment.writeArtifacts(result, summary, fig, cfg);
            second = uav.experiment.writeArtifacts(result, summary, fig, cfg);

            testCase.verifyEqual(first, second);
            paths = string(struct2cell(first));
            testCase.verifyTrue(all(isfile(paths)));
            testCase.verifyTrue(endsWith(first.mat, '.mat'));
            testCase.verifyTrue(endsWith(first.csv, '.csv'));
            testCase.verifyTrue(endsWith(first.fig, '.fig'));
            testCase.verifyTrue(endsWith(first.pdf, '.pdf'));
            testCase.verifyTrue(endsWith(first.png, '.png'));
            testCase.verifyTrue(endsWith(first.log, '.json'));
            loaded = load(first.mat, 'result');
            requiredMetadata = {'gitCommit', 'matlabVersion', 'toolboxes', ...
                'profile', 'experimentName', 'baseSeed', 'generatedAt', ...
                'codeVersion', 'actualCounts', 'seedScheme'};
            testCase.verifyTrue(all(isfield(loaded.result.metadata, requiredMetadata)));
            testCase.verifyTrue(all(isfield(loaded.result.metadata.toolboxes, ...
                {'Name', 'Version', 'Release', 'Date'})));
            testCase.verifyEqual(loaded.result.cfg, cfg);

            existingFiles = sort(string({dir(fullfile(cfg.output.rootDir, '**', '*')).name}));
            noOverwriteCfg = cfg;
            noOverwriteCfg.output.overwrite = false;
            testCase.verifyError(@() uav.experiment.writeArtifacts( ...
                result, summary, fig, noOverwriteCfg), ...
                'uav:experiment:ArtifactExists');
            testCase.verifyEqual(sort(string({dir(fullfile( ...
                cfg.output.rootDir, '**', '*')).name})), existingFiles);

            failingFigure = figure('Visible', 'off');
            failingAxes = axes('Parent', failingFigure);
            plot(failingAxes, [0, 1], [2, 1]);
            testCase.verifyError(@() uav.experiment.finalize(result, summary, ...
                failingFigure, noOverwriteCfg, struct('trials', 1), tic), ...
                'uav:experiment:ArtifactExists');
            testCase.verifyFalse(isgraphics(failingFigure));
        end

        function failedBundlePublicationRestoresPreviousFiles(testCase)
            root = tempname;
            mkdir(root);
            testCase.addTeardown(@() cleanupDir(root));
            staged = struct( ...
                'a', string(fullfile(root, 'stage-a.txt')), ...
                'b', string(fullfile(root, 'stage-b.txt')));
            artifacts = struct( ...
                'a', string(fullfile(root, 'final-a.txt')), ...
                'b', string(fullfile(root, 'missing', 'final-b.txt')));
            backups = struct( ...
                'a', string(fullfile(root, 'backup-a.txt')), ...
                'b', string(fullfile(root, 'backup-b.txt')));
            writeText(staged.a, 'new-a');
            writeText(staged.b, 'new-b');
            writeText(artifacts.a, 'old-a');

            testCase.verifyError(@() uav.experiment.publishStagedBundle( ...
                staged, artifacts, backups), ...
                'uav:experiment:ArtifactPublishFailure');

            testCase.verifyEqual(string(fileread(artifacts.a)), "old-a");
            testCase.verifyFalse(isfile(artifacts.b));
            testCase.verifyFalse(any(isfile(string(struct2cell(staged)))));
            testCase.verifyFalse(any(isfile(string(struct2cell(backups)))));
        end

        function postCommitBackupCleanupFailureKeepsNewBundle(testCase)
            root = tempname;
            mkdir(root);
            testCase.addTeardown(@() cleanupDir(root));
            staged = struct( ...
                'a', string(fullfile(root, 'stage-a.txt')), ...
                'b', string(fullfile(root, 'stage-b.txt')));
            artifacts = struct( ...
                'a', string(fullfile(root, 'final-a.txt')), ...
                'b', string(fullfile(root, 'final-b.txt')));
            backups = struct( ...
                'a', string(fullfile(root, 'backup-a.txt')), ...
                'b', string(fullfile(root, 'backup-b.txt')));
            writeText(staged.a, 'new-a');
            writeText(staged.b, 'new-b');
            writeText(artifacts.a, 'old-a');
            writeText(artifacts.b, 'old-b');
            operations = struct( ...
                'move', @moveTestFile, ...
                'delete', @(pathValue) deleteWithFailure( ...
                pathValue, backups.b));

            testCase.verifyWarning(@() uav.experiment.publishStagedBundle( ...
                staged, artifacts, backups, operations), ...
                'uav:experiment:BackupCleanupFailure');

            testCase.verifyEqual(string(fileread(artifacts.a)), "new-a");
            testCase.verifyEqual(string(fileread(artifacts.b)), "new-b");
            testCase.verifyFalse(isfile(backups.a));
            testCase.verifyEqual(string(fileread(backups.b)), "old-b");
            testCase.verifyFalse(any(isfile(string(struct2cell(staged)))));
        end

        function failedRollbackPreservesRecoverableBackup(testCase)
            root = tempname;
            mkdir(root);
            testCase.addTeardown(@() cleanupDir(root));
            staged = struct( ...
                'a', string(fullfile(root, 'stage-a.txt')), ...
                'b', string(fullfile(root, 'stage-b.txt')));
            artifacts = struct( ...
                'a', string(fullfile(root, 'final-a.txt')), ...
                'b', string(fullfile(root, 'missing', 'final-b.txt')));
            backups = struct( ...
                'a', string(fullfile(root, 'backup-a.txt')), ...
                'b', string(fullfile(root, 'backup-b.txt')));
            writeText(staged.a, 'new-a');
            writeText(staged.b, 'new-b');
            writeText(artifacts.a, 'old-a');
            operations = struct( ...
                'move', @(source, destination) moveWithRestoreFailure( ...
                source, destination, backups.a), ...
                'delete', @deleteTestFile);

            testCase.verifyError(@() uav.experiment.publishStagedBundle( ...
                staged, artifacts, backups, operations), ...
                'uav:experiment:ArtifactRollbackFailure');

            testCase.verifyFalse(isfile(artifacts.a));
            testCase.verifyEqual(string(fileread(backups.a)), "old-a");
            testCase.verifyFalse(any(isfile(string(struct2cell(staged)))));
        end

        function fiveSmokeEntrypointsProduceEveryArtifactType(testCase)
            cfg = compactConfig(testCase);
            configuredRoot = string(cfg.output.rootDir);
            entryPoints = uav.experiment.entryPoints();
            expected = ["exp01_model_validation", "exp02_jitter_sweep", ...
                "exp03_snapshot_sweep", "exp04_dynamic_stability", "exp05_ablation"];
            testCase.verifyEqual(entryPoints, expected);

            for index = 1:numel(entryPoints)
                result = feval(entryPoints(index), cfg);
                testCase.verifyEqual(result.metadata.experimentName, entryPoints(index));
                testCase.verifyGreaterThan(result.cfg.beamforming.epsilonS, 0);
                testCase.verifyEqual(result.calibration.scope, ...
                    "fixed-once-per-experiment-configuration");
                testCase.verifyEqual(result.cfg.beamforming.epsilonS, ...
                    result.calibration.epsilonS, 'AbsTol', 0);
                artifactPaths = string(struct2cell(rmfield(result.artifacts, 'status')));
                testCase.verifyTrue(all(isfile(artifactPaths)));
                testCase.verifyTrue(all(startsWith(artifactPaths, configuredRoot)));
                testCase.verifyEqual(result.artifacts.status, "complete");
                loaded = load(result.artifacts.mat, 'result');
                testCase.verifyTrue(isfield(loaded.result, 'raw'));
                testCase.verifyTrue(isfield(loaded.result, 'summary'));
                if isfield(loaded.result.raw, 'sinrLinear')
                    testCase.verifyTrue(all(isfinite(loaded.result.raw.sinrLinear), 'all'));
                end
                switch entryPoints(index)
                    case "exp01_model_validation"
                        testCase.verifyTrue(isfield(loaded.result.summary, ...
                            'zeroJitterMaximumSteeringDifference'));
                    case "exp02_jitter_sweep"
                        testCase.verifyTrue(all(isfield(loaded.result.summary, ...
                            {'sinr', 'desiredGain', 'interferenceLeakage'})));
                    case "exp03_snapshot_sweep"
                        testCase.verifyTrue(all(isfield(loaded.result.summary, ...
                            {'incmNmse', 'incmMinimumEigenvalue', ...
                            'incmConditionNumber', 'runtime'})));
                    case "exp04_dynamic_stability"
                        testCase.verifyTrue(all(isfield(loaded.result.summary, ...
                            {'outageProbability', 'periodMedianDesiredGain', ...
                            'periodMedianInterferenceLeakage', 'runtime'})));
                    case "exp05_ablation"
                        testCase.verifyTrue(isfield(loaded.result, 'switchDiagnostics'));
                end
            end
        end

        function entrypointsContainNoRuntimeNumericAssignments(testCase)
            rootDir = fileparts(fileparts(mfilename('fullpath')));
            names = uav.experiment.entryPoints();
            for index = 1:numel(names)
                source = fileread(fullfile(rootDir, 'experiments', names(index) + ".m"));
                testCase.verifyEmpty(regexp(source, '=\s*[-+]?\d', 'once'));
            end
        end

        function epsilonCalibrationOccursOnlyInExperimentPreparation(testCase)
            rootDir = fileparts(fileparts(mfilename('fullpath')));
            packageDir = fullfile(rootDir, 'src', '+uav', '+experiment');
            sources = dir(fullfile(packageDir, '*.m'));
            occurrenceCount = 0;
            occurrenceFiles = strings(0, 1);
            for index = 1:numel(sources)
                source = string(fileread(fullfile(packageDir, sources(index).name)));
                currentCount = count(source, "calibrateEpsilon(");
                occurrenceCount = occurrenceCount + currentCount;
                if currentCount > 0
                    occurrenceFiles(end + 1, 1) = string(sources(index).name); %#ok<AGROW>
                end
            end
            testCase.verifyEqual(occurrenceCount, 1);
            testCase.verifyEqual(occurrenceFiles, "prepare.m");

            runners = dir(fullfile(packageDir, 'runExp*.m'));
            testCase.verifyEqual(numel(runners), 5);
            for index = 1:numel(runners)
                source = string(fileread(fullfile(packageDir, runners(index).name)));
                testCase.verifyEqual(count(source, ...
                    "uav.experiment.prepare("), 1);
            end
        end

        function generatedArtifactsAreIgnoredByGit(testCase)
            rootDir = fileparts(fileparts(mfilename('fullpath')));
            candidates = { ...
                'results/raw/task5-ignore-check.mat', ...
                'results/tables/task5-ignore-check.csv', ...
                'results/figures/task5-ignore-check.fig', ...
                'results/logs/task5-ignore-check.json'};
            for index = 1:numel(candidates)
                [status, ~] = system(sprintf('git -C "%s" check-ignore -q "%s"', ...
                    rootDir, candidates{index}));
                testCase.verifyEqual(status, 0);
            end
        end
    end
end

function cfg = compactConfig(testCase)
cfg = config('smoke');
cfg.timeline.K = 2;
cfg.monteCarlo.nTrials = 1;
cfg.spectrum.gridSize = 9;
cfg.beamforming.calibration.sampleCount = 16;
cfg.experiments.jitterSweepDeg = [0, 1];
cfg.experiments.snapshotSweep = [2, 4];
cfg.experiments.periodCounts = struct( ...
    'modelValidation', 2, ...
    'jitterSweep', 1, ...
    'snapshotSweep', 1, ...
    'dynamicStability', 2, ...
    'ablation', 1);
cfg.experiments.bootstrapReplicates = 12;
cfg = uav.config.validate(cfg);
outputRoot = tempname;
mkdir(outputRoot);
testCase.addTeardown(@() cleanupDir(outputRoot));
cfg.output.rootDir = outputRoot;
cfg.output.rawDir = fullfile(outputRoot, 'raw');
cfg.output.tableDir = fullfile(outputRoot, 'tables');
cfg.output.figureDir = fullfile(outputRoot, 'figures');
cfg.output.logDir = fullfile(outputRoot, 'logs');
cfg.output.cacheDir = fullfile(outputRoot, 'cache');
cfg.output.overwrite = true;
end

function closeIfValid(fig)
if isgraphics(fig)
    close(fig);
end
end

function cleanupDir(pathValue)
if isfolder(pathValue)
    rmdir(pathValue, 's');
end
end

function writeText(pathValue, value)
fileId = fopen(pathValue, 'w');
if fileId < 0
    error('uav:test:FixtureWriteFailure', 'Could not create test fixture.');
end
cleanup = onCleanup(@() fclose(fileId));
fwrite(fileId, value, 'char');
end

function moveTestFile(source, destination)
[moved, message] = movefile(source, destination, 'f');
if ~moved
    error('uav:experiment:ArtifactPublishFailure', ...
        'Could not publish %s: %s', destination, message);
end
end

function moveWithRestoreFailure(source, destination, failingSource)
if string(source) == string(failingSource)
    error('uav:test:InjectedRestoreFailure', 'Injected restore failure.');
end
moveTestFile(source, destination);
end

function deleteWithFailure(pathValue, failingPath)
if string(pathValue) == string(failingPath)
    error('uav:test:InjectedCleanupFailure', 'Injected cleanup failure.');
end
deleteTestFile(pathValue);
end

function deleteTestFile(pathValue)
delete(pathValue);
end
