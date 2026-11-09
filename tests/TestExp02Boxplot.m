classdef TestExp02Boxplot < matlab.unittest.TestCase
    % 实验2箱线图的统计口径、退化场景及原入口导出验证。
    methods (TestMethodSetup)
        function addRuntimePaths(testCase)
            rootDir = fileparts(fileparts(mfilename('fullpath')));
            for pathValue = {rootDir, fullfile(rootDir, 'src'), fullfile(rootDir, 'experiments')}
                testCase.applyFixture(matlab.unittest.fixtures.PathFixture(pathValue{1}));
            end
        end
    end
    methods (Test)
        function retainsTrialUnitAndAllOutliers(testCase)
            cfg = config("smoke");
            cfg.plot.visible = false;
            % 本用例专门验证自动轴限，不依赖用户为正式实验手动设置的范围。
            cfg.plot.exp02Box.yLimits = [];
            % 每行一个试次，两周期相差0.02；最后一个试次为离群样本。
            centers = [0.10; 0.11; 0.12; 0.13; 0.14; 0.15; 0.16; 1.20];
            values = [centers - 0.01, centers + 0.01];
            raw = cat(3, zeros(size(values)), values);
            [fig, stats] = uav.plot.jitterMismatchBoxplot([0, 2], raw, cfg);
            cleanup = onCleanup(@() close(fig)); %#ok<NASGU>
            testCase.verifyEqual(stats.trialMedian, [zeros(8, 1), centers], 'AbsTol', 1e-14);
            testCase.verifyEqual(stats.quartiles, ...
                quantile([zeros(8, 1), centers], [0.25, 0.5, 0.75], 1), ...
                'AbsTol', 1e-14);
            testCase.verifyEqual(stats.whiskers(:, 2), [0.10; 0.16], 'AbsTol', 1e-14);
            testCase.verifyEqual(stats.outliers{2}, 1.20, 'AbsTol', 1e-14);
            limits = ylim(fig.CurrentAxes);
            testCase.verifyGreaterThanOrEqual(limits(2), max(centers));
            testCase.verifyEqual(stats.trialCount, 8);
            testCase.verifyEqual(stats.periodCount, 2);
            testCase.verifyNumElements(findobj(fig, 'Tag', 'MismatchBox'), 2);
            testCase.verifyEmpty(findobj(fig, 'Type', 'errorbar'));
            chartLegend = findobj(fig, 'Type', 'legend');
            testCase.verifyNumElements(chartLegend, 1);
            testCase.verifyEqual(chartLegend.FontSize, cfg.plot.fontSize);
            testCase.verifyEqual(string(chartLegend.String), ...
                "箱体：25%–75%分位");
        end

        function singleTrialZeroMismatchIsValid(testCase)
            cfg = config("smoke");
            cfg.plot.visible = false;
            cfg.plot.exp02Box.whiskerIqrMultiplier = 2;
            [fig, stats] = uav.plot.jitterMismatchBoxplot(0, 0, cfg);
            cleanup = onCleanup(@() close(fig)); %#ok<NASGU>
            testCase.verifyEqual(stats.quartiles, zeros(3, 1));
            testCase.verifyEqual(stats.whiskers, zeros(2, 1));
            testCase.verifyEmpty(stats.outliers{1});
            chartLegend = findobj(fig, 'Type', 'legend');
            testCase.verifyEqual(string(chartLegend.String), "箱体：25%–75%分位");
            testCase.verifyGreaterThan(diff(xlim(fig.CurrentAxes)), 0);
            testCase.verifyError(@() uav.plot.jitterMismatchBoxplot( ...
                0, NaN, cfg), 'uav:plot:InvalidMismatchSamples');
        end

        function experimentTwoExportsBoxplotAndKeepsStatistics(testCase)
            % 用最小的两扫描点实验确认新绘图接入原结果导出链路。
            cfg = config("smoke");
            cfg.monteCarlo.nTrials = 2;
            cfg.spectrum.gridSize = 9;
            cfg.beamforming.calibration.sampleCount = 16;
            cfg.experiments.jitterSweepDeg = [0, 1];
            cfg.experiments.periodCounts.jitterSweep = 2;
            cfg.experiments.exp02MethodNames = "oracle";
            cfg.experiments.bootstrapReplicates = 12;
            cfg.plot.visible = false;
            cfg.plot.closeAfterExport = true;
            outputRoot = tempname;
            mkdir(outputRoot);
            cleanupDir = onCleanup(@() rmdir(outputRoot, 's')); %#ok<NASGU>
            cfg.output.rootDir = outputRoot;
            cfg.output.rawDir = fullfile(outputRoot, 'raw');
            cfg.output.tableDir = fullfile(outputRoot, 'tables');
            cfg.output.figureDir = fullfile(outputRoot, 'figures');
            cfg.output.logDir = fullfile(outputRoot, 'logs');
            cfg.output.cacheDir = fullfile(outputRoot, 'cache');
            result = exp02_jitter_sweep(cfg);
            box = result.summary.desiredSteeringMismatch.boxplot;
            expected = reshape(median( ...
                result.raw.desiredSteeringMismatchNormalized, 2), 2, 2);
            testCase.verifyEqual(box.trialMedian, expected);
            testCase.verifyEqual(box.quartiles(2, :).', ...
                result.summary.desiredSteeringMismatch.median, 'AbsTol', 1e-14);
            paths = string(struct2cell(rmfield(result.artifacts, 'status')));
            testCase.verifyTrue(all(isfile(paths)));
            [~, baseName] = fileparts(result.artifacts.png);
            testCase.verifyEqual(string(baseName), "exp02_jitter_sweep_smoke");
        end
    end
end
