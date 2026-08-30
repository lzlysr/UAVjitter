classdef TestTask7GuLeshem2012Rec < matlab.unittest.TestCase
    %TESTTASK7GULESHEM2012REC Gu–Leshem 2012 Rec文献复现及UPA迁移测试。

    methods (TestMethodSetup)
        function addRuntimePaths(testCase)
            % 测试直接调用工程包函数和公开实验入口。
            rootDir = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(rootDir));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(rootDir, 'src')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(rootDir, 'experiments')));
        end
    end

    methods (Test)
        % 若实现错误地加入加载、噪声单位阵或PSD投影，本手算矩阵将不相等。
        function reconstructionMatchesPaperEquationsSixAndSeven(testCase)
            sampleCovariance = diag([2, 3]);
            scanManifold = [1, 1; 1, -1];
            quadratureWeights = [0.2, 0.3];
            expectedSpectrum = [1.2, 1.2];
            expectedIncm = [0.6, -0.12; -0.12, 0.6];

            [actualIncm, diagnostics] = uav.incm.guLeshem2012Rec( ...
                sampleCovariance, scanManifold, quadratureWeights);

            testCase.verifyEqual(actualIncm, expectedIncm, 'AbsTol', 1e-12);
            testCase.verifyEqual(diagnostics.caponPower, ...
                expectedSpectrum, 'AbsTol', 1e-12);
            testCase.verifyEqual(diagnostics.loading, 0, 'AbsTol', 0);
            testCase.verifyFalse(diagnostics.usedPseudoInverse);
        end

        % 场景测试锁定论文例1的阵列、角域、方向失配和积分总测度。
        function paperExampleOneScenarioUsesReportedParameters(testCase)
            cfg = config('smoke');
            paper = cfg.literature.guLeshem2012;
            stream = uav.util.makeStream(cfg.random.baseSeed, 710001);

            obs = uav.scenario.generateGuLeshem2012Trial( ...
                cfg, paper.fixedSnrDb, paper.fixedSnapshotCount, stream);

            testCase.verifyEqual(paper.sensorCount, 10);
            testCase.verifyEqual(paper.spacingOverLambda, 0.5);
            testCase.verifyEqual(paper.nominalDesiredDeg, 5);
            testCase.verifyEqual(paper.nominalInterfererDeg, [-50, -20]);
            testCase.verifyEqual(paper.interfererInrDb, [30, 30]);
            testCase.verifyEqual(paper.desiredSectorDeg, [0, 10]);
            testCase.verifyEqual(paper.lookMismatchHalfWidthDeg, 4);
            testCase.verifySize(obs.X, [10, paper.fixedSnapshotCount]);
            testCase.verifySize(obs.Rx, [10, 10]);
            testCase.verifyGreaterThanOrEqual(obs.actualAnglesDeg, ...
                [1, -54, -24]);
            testCase.verifyLessThanOrEqual(obs.actualAnglesDeg, ...
                [9, -46, -16]);
            testCase.verifyFalse(any(obs.scanAnglesDeg > 0 & ...
                obs.scanAnglesDeg < 10));
            testCase.verifyEqual(sum(obs.scanWeights), ...
                deg2rad(170), 'AbsTol', 1e-12);
            testCase.verifyEqual(vecnorm(obs.truth.trueASVs), ...
                sqrt(10) * ones(1, 3), 'AbsTol', 1e-12);
        end

        % 比较器只复用原论文实际比较过且工程已有的ESP与WCO。
        function paperRunnerComparesOnlyOverlappingImplementedMethods(testCase)
            cfg = config('smoke');
            paper = cfg.literature.guLeshem2012;
            obs = uav.scenario.generateGuLeshem2012Trial(cfg, ...
                paper.fixedSnrDb, paper.fixedSnapshotCount, ...
                uav.util.makeStream(cfg.random.baseSeed, 720001));

            record = uav.experiment.runGuLeshem2012Methods( ...
                obs, cfg, paper.methodNames);

            testCase.verifyEqual(record.methodNames, ...
                ["oracle", "esp", "wco", "gu_leshem2012_rec"]);
            testCase.verifyTrue(all(isfinite(record.sinrLinear)));
            testCase.verifyGreaterThanOrEqual(record.sinrLinear(1), ...
                max(record.sinrLinear(2:end)) - 1e-10);
            testCase.verifyEqual(record.diagnostics{1}.truthAccess, ...
                "oracle-only");
            testCase.verifyEqual( ...
                record.diagnostics{4}.reconstruction.loading, 0, 'AbsTol', 0);
            testCase.verifyFalse( ...
                record.diagnostics{4}.reconstruction.usedPseudoInverse);
            testCase.verifyEqual(record.snapshotSignatures, ...
                repmat(record.snapshotSignatures(1, :), 4, 1), 'AbsTol', 0);
        end

        % 独立复现实验必须同时输出论文图1(a)与图1(c)对应的两组扫描。
        function paperReproductionExperimentPublishesBothSweeps(testCase)
            cfg = compactLiteratureConfig(testCase);
            cfg.literature.guLeshem2012.snrSweepDb = [-30, 20];
            cfg.literature.guLeshem2012.snapshotSweep = [20, 30];

            result = exp07_gu_leshem2012_reproduction(cfg);

            testCase.verifyEqual(result.methodNames, ...
                cfg.literature.guLeshem2012.methodNames);
            testCase.verifySize(result.raw.snrSinrLinear, [1, 2, 4]);
            testCase.verifySize(result.raw.snapshotSinrLinear, [1, 2, 4]);
            testCase.verifySize(result.summary.snr.meanDb, [2, 4]);
            testCase.verifySize(result.summary.snapshot.meanDb, [2, 4]);
            testCase.verifyFalse(any(uav.experiment.entryPoints() == ...
                "exp07_gu_leshem2012_reproduction"));
            artifactFields = {'mat', 'csv', 'fig', 'pdf', 'png', 'log'};
            artifactPaths = strings(numel(artifactFields), 1);
            for artifactIndex = 1:numel(artifactFields)
                artifactPaths(artifactIndex) = ...
                    result.artifacts.(artifactFields{artifactIndex});
            end
            testCase.verifyTrue(all(isfile(artifactPaths)));
        end

        % UPA迁移只改变阵列与二维求积，不加入论文之外的动态或时序算法。
        function upaTransferReportsGridAndSectorDiagnostics(testCase)
            cfg = compactLiteratureConfig(testCase);
            cfg.literature.guLeshem2012.upaTransferPeriods = 1;
            cfg.literature.guLeshem2012.upaTransferJitterDeg = [0, 8];
            cfg.literature.guLeshem2012.upaTransferGridSizes = [9, 11];

            result = exp08_gu_leshem2012_upa_transfer(cfg);

            testCase.verifySize(result.raw.sinrLinear, [1, 1, 2, 2, 4]);
            testCase.verifySize(result.summary.sinr.meanDb, [2, 2, 4]);
            testCase.verifySize(result.raw.incmNmse, [1, 1, 2, 2]);
            testCase.verifySize(result.raw.normalizedEigenvalueProfileError, ...
                [1, 1, 2, 2]);
            testCase.verifySize(result.raw.estimatedConditionNumber, ...
                [1, 1, 2, 2]);
            testCase.verifyTrue(all(isfinite( ...
                result.summary.normalizedEigenvalueProfileErrorMean), 'all'));
            testCase.verifyTrue(all( ...
                result.summary.estimatedConditionNumberMedian >= 1, 'all'));
            testCase.verifySize(result.raw.desiredOutsideExclusion, ...
                [1, 1, 2, 2]);
            testCase.verifySize(result.raw.interfererInsideExclusion, ...
                [1, 1, 2, 2, 2]);
            testCase.verifyFalse(any(uav.experiment.entryPoints() == ...
                "exp08_gu_leshem2012_upa_transfer"));
            artifactFields = {'mat', 'csv', 'fig', 'pdf', 'png', 'log'};
            for artifactIndex = 1:numel(artifactFields)
                testCase.verifyTrue(isfile( ...
                    result.artifacts.(artifactFields{artifactIndex})));
            end
        end

        % 原式要求SCM可逆，且文献比较方法不得混入论文未比较的工程方法。
        function literatureConfigurationRejectsOutOfScopeOrSingularCases(testCase)
            cfg = config('smoke');

            invalidMethod = cfg;
            invalidMethod.literature.guLeshem2012.methodNames = ...
                ["oracle", "dl_mvdr", "gu_leshem2012_rec"];
            testCase.verifyError(@() uav.config.validate(invalidMethod), ...
                'uav:config:InvalidExperimentParameter');

            missingReconstruction = cfg;
            missingReconstruction.literature.guLeshem2012.methodNames = ...
                ["oracle", "esp", "wco"];
            testCase.verifyError(@() uav.config.validate(missingReconstruction), ...
                'uav:config:InvalidExperimentParameter');

            % 用户可能在config()返回后直接修改字段；实验入口也必须自校验。
            directInvalid = compactLiteratureConfig(testCase);
            directInvalid.literature.guLeshem2012.methodNames = ...
                ["oracle", "esp", "wco"];
            testCase.verifyError( ...
                @() exp07_gu_leshem2012_reproduction(directInvalid), ...
                'uav:config:InvalidExperimentParameter');

            singularSweep = cfg;
            singularSweep.literature.guLeshem2012.snapshotSweep = 9;
            testCase.verifyError(@() uav.config.validate(singularSweep), ...
                'uav:config:InvalidExperimentParameter');

            unsupportedGrid = cfg;
            unsupportedGrid.literature.guLeshem2012.upaTransferGridSizes = 10;
            testCase.verifyError(@() uav.config.validate(unsupportedGrid), ...
                'uav:config:InvalidExperimentParameter');
        end
    end
end

function cfg = compactLiteratureConfig(testCase)
%COMPACTLITERATURECONFIG 构造文献复现实验的单试次临时输出配置。
cfg = config('smoke');
cfg.literature.guLeshem2012.monteCarloTrials = 1;
cfg.experiments.bootstrapReplicates = 8;
cfg.plot.visible = false;
cfg.plot.closeAfterExport = true;
outputRoot = tempname;
mkdir(outputRoot);
testCase.addTeardown(@() cleanupLiteratureDirectory(outputRoot));
cfg.output.rootDir = outputRoot;
cfg.output.rawDir = fullfile(outputRoot, 'raw');
cfg.output.tableDir = fullfile(outputRoot, 'tables');
cfg.output.figureDir = fullfile(outputRoot, 'figures');
cfg.output.logDir = fullfile(outputRoot, 'logs');
cfg.output.cacheDir = fullfile(outputRoot, 'cache');
end

function cleanupLiteratureDirectory(pathValue)
%CLEANUPLITERATUREDIRECTORY 删除测试专用且可重建的结果目录。
if isfolder(pathValue)
    rmdir(pathValue, 's');
end
end
