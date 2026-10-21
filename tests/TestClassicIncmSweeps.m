classdef TestClassicIncmSweeps < matlab.unittest.TestCase
    %TESTCLASSICINCMSWEEPS 验证原static复用及两张经典方法扫描图的数据契约。
    methods (TestMethodSetup)
        function preparePaths(testCase)
            root = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root, 'src')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root, 'experiments')));
        end
    end
    methods (Test)
        function staticUsesExistingReconstructionAndStandardMvdr(testCase)
            cfg = config("smoke");
            scene = uav.scenario.generate(cfg, uav.util.makeStream(cfg.random.baseSeed, 1));
            obs = uav.scenario.generatePeriod(cfg, scene, 1, ...
                uav.util.makeStream(cfg.random.baseSeed, 2));
            safe = uav.experiment.algorithmObservation(obs);
            [matrix, ~] = uav.incm.reconstruct("static", safe, [], cfg);
            expected = uav.beamforming.solve("mvdr", matrix, ...
                cfg.derived.nominalSteeringReference, cfg);
            record = uav.experiment.runClassicMethods(obs, [], cfg, ...
                ["oracle", "static", "CovFit"]);
            testCase.verifyEqual(record.weights(:, 2), expected, 'AbsTol', 1e-12);
            testCase.verifyEqual(record.diagnostics{2}.beamforming.method, "mvdr");
            testCase.verifyTrue(all(isfinite(record.sinrLinear)));
            testCase.verifyEqual(record.snapshotSignatures(1, :), record.snapshotSignatures(2, :));
            testCase.verifyEqual(record.snapshotSignatures(1, :), record.snapshotSignatures(3, :));
        end
        function snapshotSelectionRejectsDuplicatesAndUnknownNames(testCase)
            cfg = config("smoke");
            cfg.experiments.classicSnapshotMethodNames = ["static", "oracle"];
            checked = uav.config.validate(cfg);
            testCase.verifyEqual(checked.experiments.classicSnapshotMethodNames, ["static", "oracle"]);
            cfg.experiments.classicSnapshotMethodNames = ["static", "static"];
            testCase.verifyError(@() uav.config.validate(cfg), 'uav:config:InvalidExperimentParameter');
            cfg.experiments.classicSnapshotMethodNames = "unknown";
            testCase.verifyError(@() uav.config.validate(cfg), 'uav:config:InvalidExperimentParameter');
        end
        function bothSweepsExportPairedSelectedMethodResults(testCase)
            cfg = config("smoke");
            folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            cfg.output.rawDir = fullfile(folder.Folder, 'raw');
            cfg.output.figureDir = fullfile(folder.Folder, 'figures');
            cfg.output.tableDir = fullfile(folder.Folder, 'tables');
            cfg.output.logDir = fullfile(folder.Folder, 'logs');
            % 人工小夹具只验证执行契约，不用于汇报性能结论。
            cfg.monteCarlo.nTrials = 2;
            cfg.experiments.jitterSweepDeg = [0, 2];
            cfg.experiments.snapshotSweep = [8, 32];
            cfg.experiments.bootstrapReplicates = 10;
            cfg.plot.visible = false;
            cfg.plot.closeAfterExport = true;
            originalStd = cfg.attitude.stdDeg;
            jitter = exp06_classic_method_jitter_sweep(cfg);
            snapshots = exp07_classic_method_snapshot_sweep(cfg);
            for result = {jitter, snapshots}
                value = result{1};
                testCase.verifyEqual(value.status, "complete");
                testCase.verifyEqual(size(value.raw.sinrLinear), [2, 2, 2, numel(value.methodNames)]);
                testCase.verifyTrue(all(isfinite(value.raw.sinrLinear), 'all'));
                testCase.verifyTrue(isfile(value.artifacts.mat));
                testCase.verifyTrue(isfile(value.artifacts.fig));
                testCase.verifyTrue(isfile(value.artifacts.pdf));
                testCase.verifyTrue(isfile(value.artifacts.png));
                testCase.verifyTrue(isfile(value.artifacts.csv));
                for method = 2:numel(value.methodNames)
                    testCase.verifyEqual(value.raw.snapshotSignatures(:, :, :, 1, :), ...
                        value.raw.snapshotSignatures(:, :, :, method, :));
                end
            end
            testCase.verifyEqual(snapshots.cfg.attitude.stdDeg, originalStd);
            testCase.verifyEqual(snapshots.scanPoints, [8, 32]);
            testCase.verifyNotEqual(jitter.seedScheme{1}.scenarioSubstreams, ...
                snapshots.seedScheme{1}.scenarioSubstreams);
            testCase.verifyEqual(jitter.methodNames, snapshots.methodNames);
        end
    end
end
