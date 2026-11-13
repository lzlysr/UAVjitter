classdef TestClassicCalibrationGate < matlab.unittest.TestCase
    % 检查经典扫描按所选方法标定半径，不让无关WCO阻断标准MVDR对照。
    methods (TestMethodSetup)
        function addProjectPaths(testCase)
            root = fileparts(fileparts(mfilename('fullpath')));
            for path = {root, fullfile(root, 'src'), fullfile(root, 'experiments')}
                testCase.applyFixture(matlab.unittest.fixtures.PathFixture(path{1}));
            end
        end
    end
    methods (Test)
        function withoutWcoBothArraysCompleteAllThreeScans(testCase)
            folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            for shape = {[9, 9], [9, 1]}
                cfg = compactConfig(shape{1}, folder.Folder);
                % 9×9下原设计半径确实不可行，跳过此标定才是本测试的关键。
                if cfg.derived.arrayType == "UPA"
                    testCase.verifyError(@() uav.beamforming.calibrateEpsilon(cfg), ...
                        'uav:beamforming:InfeasibleCalibration');
                end
                for kind = ["jitter", "snapshots", "snr"]
                    result = uav.experiment.runClassicSweep(cfg, kind);
                    testCase.verifyEqual(result.status, "complete");
                    testCase.verifyEqual(result.cfg.timeline.L, 96);
                    testCase.verifyEqual(result.cfg.beamforming.epsilonS, ...
                        cfg.beamforming.epsilonS);
                    testCase.verifyEqual(result.calibration.epsilonS, ...
                        cfg.beamforming.epsilonS);
                    testCase.verifyFalse(isfield(result.calibration, 'diagnostics'));
                    testCase.verifyTrue(contains(result.calibration.scope, "未选择WCO"));
                    testCase.verifyTrue(all(result.raw.sinrLinear > 0, 'all'));
                    testCase.verifyTrue(all(isfinite(result.raw.sinrLinear), 'all'));
                    testCase.verifyTrue(all(result.raw.methodStatus == "ok", 'all'));
                    testCase.verifyTrue(isfile(result.artifacts.mat));
                end
            end
        end
        function selectingWcoStillRejectsInfeasibleUpaRadius(testCase)
            folder = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            cfg = compactConfig([9, 9], folder.Folder);
            cfg.experiments.classicJitterMethodNames = ["oracle", "wco"];
            cfg.experiments.classicSnapshotMethodNames = ["oracle", "wco"];
            cfg.experiments.classicSnrMethodNames = ["oracle", "wco"];
            for kind = ["jitter", "snapshots", "snr"]
                testCase.verifyError(@() uav.experiment.runClassicSweep(cfg, kind), ...
                    'uav:beamforming:InfeasibleCalibration');
            end
        end
    end
end

function cfg = compactConfig(shape, outputRoot)
% 人工单点夹具仅用于流程回归，不据此评价算法性能，不写入正式结果目录。
cfg = config("smoke");
cfg.array.Nx = shape(1);
cfg.array.Ny = shape(2);
cfg.spectrum.gridSize = 9;
cfg.spectrum.ulaGridSize = 41;
cfg.timeline.L = 96;
cfg.beamforming.epsilonS = 0.25;
cfg.monteCarlo.nTrials = 1;
cfg.experiments.jitterSweepDeg = 0;
cfg.experiments.snapshotSweep = 96;
cfg.experiments.snrSweepDb = 10;
cfg.experiments.classicJitterMethodNames = ["oracle", "incm_integral_origin"];
cfg.experiments.classicSnapshotMethodNames = ["oracle", "incm_integral_origin"];
cfg.experiments.classicSnrMethodNames = ["oracle", "incm_integral_origin"];
cfg.experiments.periodCounts.jitterSweep = 1;
cfg.experiments.periodCounts.snapshotSweep = 1;
cfg.experiments.periodCounts.snrSweep = 1;
cfg.experiments.bootstrapReplicates = 10;
cfg.plot.visible = false;
cfg.plot.closeAfterExport = true;
for name = {'rawDir', 'figureDir', 'tableDir', 'logDir'}
    cfg.output.(name{1}) = fullfile(outputRoot, name{1});
end
cfg = uav.config.validate(cfg);
end
