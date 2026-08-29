classdef TestTask5Experiments < matlab.unittest.TestCase
    %TESTTASK5EXPERIMENTS 第五阶段：实验编排、统计、导出和故障恢复测试。
    % 这些测试锁定“方法公平比较”和“结果包完整发布”边界，不评价某个
    % 方法的科学优越性。所有产生的文件都写入临时目录并在测试后清理。
    methods (TestMethodSetup)
        % 同时加入根目录、src和experiments，以便测试公开入口和包实现。
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
        % 验证方法/消融顺序、周期数、bootstrap和绘图输出等实验配置合同。
        function configDefinesAndValidatesAllExperimentControls(testCase)
            cfg = config('smoke');

            testCase.verifyEqual(cfg.experiments.methodNames, ...
                ["oracle", "frozen", "scm", "dl", "static", "temporal"]);
            testCase.verifyFalse(isempty(cfg.experiments.exp02MethodNames));
            testCase.verifyEqual(numel(unique( ...
                cfg.experiments.exp02MethodNames)), ...
                numel(cfg.experiments.exp02MethodNames));
            testCase.verifyTrue(all(ismember( ...
                cfg.experiments.exp02MethodNames, ...
                cfg.experiments.methodNames)));
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
            invalid.experiments.exp02MethodNames = ["scm", "scm"];
            testCase.verifyError(@() uav.config.validate(invalid), ...
                'uav:config:InvalidExperimentParameter');
            invalid = cfg;
            invalid.experiments.exp02MethodNames = "unknown";
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

        % 线性域统计与共享行bootstrap必须可复现，且不能直接平均dB样本。
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

        % 试次簇重采样必须保持方法配对，四倍线性关系应在CI中精确保留。
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

        % paper的500×500数据只能按500个试次簇抽样，不能展开周期bootstrap。
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

        % 六种方法应共享一个快拍，冻结权值保持不变，并使用固定epsilonS。
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

        % 用保存的子流独立重建每个obs，逐方法核对快拍摘要以证明真实配对。
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

        % 实验采集层必须从同一姿态场景记录式（3-10）的期望导向失配，
        % 且该指标与六种方法无关，不能沿方法维重复保存。
        function collectMethodsRecordsDesiredSteeringMismatch(testCase)
            cfg = compactConfig(testCase);
            cfg.monteCarlo.nTrials = 2;
            cfg.attitude.stdDeg = [2, 1, 3];
            cfg = uav.config.validate(cfg);

            [raw, seeds] = uav.experiment.collectMethods(cfg, 2, 1, 2);

            testCase.verifySize(raw.desiredSteeringMismatchNormalized, [2, 2]);
            nominal = cfg.derived.nominalSteeringReference;
            for trialIndex = 1:cfg.monteCarlo.nTrials
                scenario = uav.scenario.generate(cfg, uav.util.makeStream( ...
                    cfg.random.baseSeed, seeds.scenarioSubstreams(trialIndex)));
                for periodIndex = 1:2
                    actual = scenario.steeringTruth(:, 1, periodIndex);
                    expected = norm(actual - nominal, 2) / norm(nominal, 2);
                    testCase.verifyEqual( ...
                        raw.desiredSteeringMismatchNormalized( ...
                        trialIndex, periodIndex), expected, 'AbsTol', 1e-14);
                end
            end
        end

        % 实验2必须把期望导向失配作为独立于方法的抖动RMS扫描指标，
        % 同时保存逐周期原始值、统计汇总和CSV列。
        function jitterSweepExportsDesiredSteeringMismatch(testCase)
            cfg = compactConfig(testCase);

            result = uav.experiment.runExp02(cfg);

            testCase.verifySize( ...
                result.raw.desiredSteeringMismatchNormalized, [1, 1, 2]);
            testCase.verifyEqual( ...
                result.raw.desiredSteeringMismatchNormalized(:, :, 1), ...
                0, 'AbsTol', 1e-14);
            testCase.verifyGreaterThan( ...
                result.raw.desiredSteeringMismatchNormalized(:, :, 2), 0);
            testCase.verifyTrue(isfield(result.summary, ...
                'desiredSteeringMismatch'));
            testCase.verifySize( ...
                result.summary.desiredSteeringMismatch.median, [2, 1]);

            exported = readtable(result.artifacts.csv);
            testCase.verifyTrue(ismember( ...
                'medianDesiredSteeringMismatchNormalized', ...
                exported.Properties.VariableNames));
            for pointIndex = 1:numel(cfg.experiments.jitterSweepDeg)
                rows = exported.jitterRmsDeg == ...
                    cfg.experiments.jitterSweepDeg(pointIndex);
                testCase.verifyEqual(unique( ...
                    exported.medianDesiredSteeringMismatchNormalized(rows)), ...
                    result.summary.desiredSteeringMismatch.median(pointIndex), ...
                    'AbsTol', 1e-14);
            end
        end

        % 非Oracle白名单应剥离truth以及未来可能增加的未批准诊断字段。
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

        % 六个消融变体必须各自只改变声明的一个配置字段。
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

        % 结果写入器应生成六种格式、完整元数据，并能确定性覆盖同名结果。
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

        % 实验2应只计算配置选择的方法，并保持用户给出的顺序；不选择
        % Oracle时也不能执行依赖“Oracle位于第一列”的上界检查。
        function experimentTwoSupportsSelectedMethodSubset(testCase)
            cfg = compactConfig(testCase);
            cfg.experiments.exp02MethodNames = ["dl", "scm"];
            % 集成测试完成导出后关闭图窗，避免批处理会话等待窗口退出。
            cfg.plot.closeAfterExport = true;
            cfg = uav.config.validate(cfg);

            [raw, ~] = uav.experiment.collectMethods(cfg, 2, 1, 1, ...
                cfg.experiments.exp02MethodNames);

            testCase.verifySize(raw.sinrLinear, [1, 1, 2]);
            testCase.verifySize(raw.snapshotSignatures, [1, 1, 2, 4]);
            testCase.verifyTrue(all(isfinite(raw.sinrLinear), 'all'));

            result = uav.experiment.runExp02(cfg);
            testCase.verifyEqual(result.methodNames, ["dl", "scm"]);
            testCase.verifySize(result.raw.sinrLinear, [1, 1, 2, 2]);
            exported = readtable(result.artifacts.csv, ...
                'TextType', 'string');
            testCase.verifyEqual(unique(exported.method, 'stable').', ...
                ["dl", "scm"]);
        end

        % 单独运行实验且未要求导出后关闭时，统一收尾模块必须保留图窗。
        function singleExperimentFinalizeKeepsFigureOpen(testCase)
            cfg = compactConfig(testCase);
            cfg.plot.closeAfterExport = false;
            fig = figure('Visible', 'off');
            plot([0, 1], [1, 2]);
            testCase.addTeardown(@() closeIfValid(fig));
            result = struct( ...
                'experimentName', "single_experiment_figure_lifecycle", ...
                'methodNames', strings(1, 0), ...
                'raw', struct('value', 1), ...
                'summary', struct('value', 1));
            summary = table(1, 'VariableNames', {'value'});

            uav.experiment.finalize(result, summary, fig, cfg, ...
                struct('trials', 1), tic);

            testCase.verifyTrue(isgraphics(fig, 'figure'));
        end

        % 发布中途失败时必须恢复旧文件并清理暂存与备份。
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

        % 到达提交点后备份清理失败不得撤销完整新结果，应保留旧备份报警。
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

        % 回滚恢复本身失败时必须保留可人工恢复的backup并报告专用错误。
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

        % 五个smoke入口必须分别生成MAT/CSV/FIG/PDF/PNG/JSON且字段完整。
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

        % 公开实验入口不得出现运行时数值赋值，防止参数绕过config。
        function entrypointsContainNoRuntimeNumericAssignments(testCase)
            rootDir = fileparts(fileparts(mfilename('fullpath')));
            names = uav.experiment.entryPoints();
            for index = 1:numel(names)
                source = fileread(fullfile(rootDir, 'experiments', names(index) + ".m"));
                testCase.verifyEmpty(regexp(source, '=\s*[-+]?\d', 'once'));
            end
        end

        % epsilon标定调用只允许出现在prepare中，五个runExp各调用prepare一次。
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

        % 生成结果和缓存必须被Git忽略，避免大型可再生文件进入版本库。
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
%COMPACTCONFIG 创建临时输出目录和小规模完整实验配置，加快集成测试。
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
%CLOSEIFVALID 清理仍存在的测试图窗。
if isgraphics(fig)
    close(fig);
end
end

function cleanupDir(pathValue)
%CLEANUPDIR 删除本测试独占的临时目录及其中可再生成文件。
if isfolder(pathValue)
    rmdir(pathValue, 's');
end
end

function writeText(pathValue, value)
%WRITETEXT 创建事务故障注入所需的短文本夹具。
fileId = fopen(pathValue, 'w');
if fileId < 0
    error('uav:test:FixtureWriteFailure', 'Could not create test fixture.');
end
cleanup = onCleanup(@() fclose(fileId));
fwrite(fileId, value, 'char');
end

function moveTestFile(source, destination)
%MOVETESTFILE 与生产move语义一致的测试操作封装。
[moved, message] = movefile(source, destination, 'f');
if ~moved
    error('uav:experiment:ArtifactPublishFailure', ...
        'Could not publish %s: %s', destination, message);
end
end

function moveWithRestoreFailure(source, destination, failingSource)
%MOVEWITHRESTOREFAILURE 在指定backup恢复时注入可控失败。
if string(source) == string(failingSource)
    error('uav:test:InjectedRestoreFailure', 'Injected restore failure.');
end
moveTestFile(source, destination);
end

function deleteWithFailure(pathValue, failingPath)
%DELETEWITHFAILURE 在指定备份清理时注入可控失败。
if string(pathValue) == string(failingPath)
    error('uav:test:InjectedCleanupFailure', 'Injected cleanup failure.');
end
deleteTestFile(pathValue);
end

function deleteTestFile(pathValue)
%DELETETESTFILE 默认测试删除操作。
delete(pathValue);
end
