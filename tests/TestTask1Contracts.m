classdef TestTask1Contracts < matlab.unittest.TestCase
    %TESTTASK1CONTRACTS 第一阶段：工程结构、配置和复现性合同测试。
    % 这些测试保证后续模块读取的是同一套参数语义，并防止profile、路径或
    % 随机流在不知情的情况下改变物理模型。
    methods (TestMethodSetup)
        % 每个测试前只加入src；PathFixture会在测试结束后自动恢复搜索路径。
        function addSourcePath(testCase)
            rootDir = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(rootDir));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(...
                fullfile(rootDir, 'src')));
        end
    end

    methods (Test)
        % 验证默认3×3阵列、两个干扰、功率和名义流形等批准的初始条件。
        function configSuppliesApprovedDefaultsAndDerivedQuantities(testCase)
            cfg = config('smoke');

            testCase.verifyEqual(cfg.array.Nx, 3);
            testCase.verifyEqual(cfg.array.Ny, 3);
            testCase.verifyEqual(cfg.array.spacingOverLambda, [0.5, 0.5]);
            testCase.verifyEqual(cfg.sources.desiredAzElDeg, [0, -35]);
            testCase.verifyEqual(cfg.sources.interfererAzElDeg, [40, -25; -55, -20]);
            testCase.verifySize(cfg.sources.interfererInrDb, [1, 2]);
            testCase.verifyTrue(all(isfinite(cfg.sources.interfererInrDb)));
            testCase.verifyEqual(cfg.signal.snrDb, 10);
            testCase.verifyEqual(cfg.timeline.dtSec, 1e-3);
            testCase.verifyGreaterThan(cfg.timeline.L, 0);
            testCase.verifyEqual(cfg.timeline.L, floor(cfg.timeline.L));
            testCase.verifyEqual(cfg.attitude.rho, [0.95, 0.95, 0.98]);
            testCase.verifyEqual(cfg.random.baseSeed, 20260827);
            testCase.verifyEqual(cfg.experiments.jitterSweepDeg(1), 0);
            testCase.verifyTrue(all(diff(cfg.experiments.jitterSweepDeg) > 0));
            testCase.verifyEqual(cfg.experiments.snapshotSweep, [2, 4, 8, 16, 32, 64]);
            testCase.verifyEqual(cfg.derived.M, 9);
            testCase.verifyEqual(cfg.derived.J, 2);
            testCase.verifyEqual(cfg.derived.noisePower, 1);
            testCase.verifyEqual(cfg.derived.desiredPower, 10, 'AbsTol', 1e-12);
            expectedInterfererPowers = cfg.signal.noisePower * ...
                10.^(cfg.sources.interfererInrDb / 10);
            testCase.verifyEqual(cfg.derived.interfererPowers, ...
                expectedInterfererPowers, ...
                'AbsTol', 1e-12);
            testCase.verifyEqual(cfg.derived.worldUnitVectors(:, 1), ...
                [0.819152044288992; 0; -0.573576436351046], 'AbsTol', 1e-12);
            testCase.verifySize(cfg.derived.nominalSteeringReference, [9, 1]);
            testCase.verifyEqual(norm(cfg.derived.nominalSteeringReference), 3, ...
                'AbsTol', 1e-12);
            testCase.verifyTrue(isfolder(cfg.output.cacheDir));
            testCase.verifyTrue(isfolder(cfg.output.rawDir));
        end

        % 验证三个profile只改变工作量/网格，物理与算法配置保持相同。
        function profilesOnlyChangeWorkloadAndGridResolution(testCase)
            smoke = config('smoke');
            preliminary = config('preliminary');
            paper = config('paper');

            testCase.verifyLessThan(smoke.monteCarlo.nTrials, ...
                preliminary.monteCarlo.nTrials);
            testCase.verifyLessThan(preliminary.monteCarlo.nTrials, ...
                paper.monteCarlo.nTrials);
            testCase.verifyLessThan(smoke.timeline.K, preliminary.timeline.K);
            testCase.verifyLessThan(preliminary.timeline.K, paper.timeline.K);
            testCase.verifyLessThan(smoke.spectrum.gridSize, preliminary.spectrum.gridSize);
            testCase.verifyLessThan(preliminary.spectrum.gridSize, paper.spectrum.gridSize);
            testCase.verifyLessThan( ...
                smoke.literature.guLeshem2012.monteCarloTrials, ...
                preliminary.literature.guLeshem2012.monteCarloTrials);
            testCase.verifyLessThan( ...
                preliminary.literature.guLeshem2012.monteCarloTrials, ...
                paper.literature.guLeshem2012.monteCarloTrials);
            testCase.verifyGreaterThan( ...
                smoke.literature.guLeshem2012.angleGridStepDeg, ...
                preliminary.literature.guLeshem2012.angleGridStepDeg);
            testCase.verifyGreaterThan( ...
                preliminary.literature.guLeshem2012.angleGridStepDeg, ...
                paper.literature.guLeshem2012.angleGridStepDeg);

            smoke = rmfield(smoke, {'profile', 'derived', 'output'});
            preliminary = rmfield(preliminary, {'profile', 'derived', 'output'});
            paper = rmfield(paper, {'profile', 'derived', 'output'});
            smoke.timeline.K = [];
            preliminary.timeline.K = [];
            paper.timeline.K = [];
            smoke.monteCarlo.nTrials = [];
            preliminary.monteCarlo.nTrials = [];
            paper.monteCarlo.nTrials = [];
            smoke.spectrum.gridSize = [];
            preliminary.spectrum.gridSize = [];
            paper.spectrum.gridSize = [];
            literatureWorkloadFields = {'monteCarloTrials', ...
                'angleGridStepDeg', 'upaTransferGridSizes', ...
                'upaTransferPeriods'};
            for index = 1:numel(literatureWorkloadFields)
                field = literatureWorkloadFields{index};
                smoke.literature.guLeshem2012.(field) = [];
                preliminary.literature.guLeshem2012.(field) = [];
                paper.literature.guLeshem2012.(field) = [];
            end
            smoke.literature.guLeshem2012.upaDiagnosis.monteCarloTrials = [];
            preliminary.literature.guLeshem2012.upaDiagnosis.monteCarloTrials = [];
            paper.literature.guLeshem2012.upaDiagnosis.monteCarloTrials = [];
            smoke.literature.guLeshem2012.upaDiagnosis.periods = [];
            preliminary.literature.guLeshem2012.upaDiagnosis.periods = [];
            paper.literature.guLeshem2012.upaDiagnosis.periods = [];
            workloadExperimentFields = {'modelValidation', 'jitterSweep', ...
                'snapshotSweep', 'dynamicStability', 'ablation'};
            for index = 1:numel(workloadExperimentFields)
                field = workloadExperimentFields{index};
                smoke.experiments.periodCounts.(field) = [];
                preliminary.experiments.periodCounts.(field) = [];
                paper.experiments.periodCounts.(field) = [];
            end
            smoke.experiments.bootstrapReplicates = [];
            preliminary.experiments.bootstrapReplicates = [];
            paper.experiments.bootstrapReplicates = [];
            testCase.verifyEqual(smoke, preliminary);
            testCase.verifyEqual(preliminary, paper);
        end

        % 未知profile必须立即报错，不能静默落到某个默认规模。
        function configRejectsUnknownProfile(testCase)
            testCase.verifyError(@() config('draft'), 'uav:config:UnknownProfile');
        end

        % 拒绝不可见方向和干扰方向数与INR长度不匹配的配置。
        function validationRejectsNonVisibleSourcesAndMismatchedPowers(testCase)
            cfg = config('smoke');
            cfg.sources.desiredAzElDeg(2) = 1;
            testCase.verifyError(@() uav.config.validate(cfg), ...
                'uav:config:SourceNotLowerHemisphere');

            cfg = config('smoke');
            cfg.sources.interfererInrDb = 20;
            testCase.verifyError(@() uav.config.validate(cfg), ...
                'uav:config:InrLengthMismatch');
        end

        % J=0是正式支持的边界场景，派生功率与流形维度必须仍然一致。
        function validationSupportsZeroInterferers(testCase)
            cfg = config('smoke');
            cfg.sources.interfererAzElDeg = zeros(0, 2);
            cfg.sources.interfererInrDb = [];

            validated = uav.config.validate(cfg);

            testCase.verifyEqual(validated.derived.J, 0);
            testCase.verifyEqual(validated.derived.interfererPowers, zeros(1, 0));
            testCase.verifySize(validated.derived.nominalSteering, [9, 1]);
            testCase.verifySize(validated.derived.nominalSteeringReference, [9, 1]);
        end

        % 对负正则项、非法阈值等算法参数进行集中范围检查。
        function validationRejectsInvalidAlgorithmParameters(testCase)
            invalidCases = { ...
                'spectrum', 'desiredExclusionHalfWidthDeg', NaN; ...
                'spectrum', 'desiredExclusionHalfWidthDeg', -1; ...
                'spectrum', 'diagonalLoading', -1; ...
                'incm', 'shrinkage', 1.01; ...
                'incm', 'eigenvalueFloorRatio', 0; ...
                'incm', 'lambdaS', -1; ...
                'incm', 'lambdaT', -Inf; ...
                'beamforming', 'diagonalLoading', NaN; ...
                'beamforming', 'epsilonS', -1};

            for index = 1:size(invalidCases, 1)
                cfg = config('smoke');
                section = invalidCases{index, 1};
                field = invalidCases{index, 2};
                cfg.(section).(field) = invalidCases{index, 3};
                testCase.verifyError(@() uav.config.validate(cfg), ...
                    'uav:config:InvalidAlgorithmParameter');
            end
        end

        % 验证Nx不等于Ny时M、源流形和参考矢量均按矩形阵列派生。
        function validationDerivesRectangularArraySteering(testCase)
            cfg = config('smoke');
            cfg.array.Nx = 2;
            cfg.array.Ny = 4;

            validated = uav.config.validate(cfg);

            testCase.verifyEqual(validated.derived.M, 8);
            testCase.verifySize(validated.derived.nominalSteering, ...
                [8, validated.derived.J + 1]);
            testCase.verifySize(validated.derived.nominalSteeringReference, [8, 1]);
        end

        % 局部随机流必须可复现，并且不得改变MATLAB全局rng序列。
        function makeStreamIsDeterministicAndDoesNotChangeGlobalRng(testCase)
            streamA = uav.util.makeStream(42, 7);
            streamB = uav.util.makeStream(42, 7);
            testCase.verifyEqual(rand(streamA, 1, 4), rand(streamB, 1, 4));

            rng(31, 'twister');
            expected = rand(1, 3);
            rng(31, 'twister');
            uav.util.makeStream(42, 8);
            testCase.verifyEqual(rand(1, 3), expected);
        end

        % startup只能添加本工程运行目录，不能执行实验或加入无关路径。
        function startupOnlyAddsRepositoryRuntimePaths(testCase)
            rootDir = fileparts(fileparts(mfilename('fullpath')));
            expected = string({fullfile(rootDir, 'src'), fullfile(rootDir, 'experiments'), ...
                fullfile(rootDir, 'tests')});
            before = string(strsplit(path, pathsep));

            startup;

            after = string(strsplit(path, pathsep));
            added = setdiff(after, before);
            testCase.verifyEmpty(setdiff(added, expected));
            testCase.verifyTrue(all(ismember(expected, after)));
        end

        % run_all允许关闭测试和实验，用于检查编排器本身而不产生副作用。
        function runAllSupportsAQuietOrchestrationModeForTests(testCase)
            result = run_all('smoke', 'RunTests', false, 'RunExperiments', false);

            testCase.verifyEqual(result.profile, "smoke");
            testCase.verifyEmpty(result.tests);
            testCase.verifyEmpty(result.experiments);
        end
    end
end
