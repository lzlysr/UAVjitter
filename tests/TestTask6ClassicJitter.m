classdef TestTask6ClassicJitter < matlab.unittest.TestCase
    %TESTTASK6CLASSICJITTER 经典波束形成方法及独立抖动扫描实验测试。

    methods (TestMethodSetup)
        function addRuntimePaths(testCase)
            % 加入根目录、源码和公开实验入口路径。
            rootDir = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(rootDir));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(rootDir, 'src')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(rootDir, 'experiments')));
        end
    end

    methods (Test)
        % 新实验必须提供完整默认集合，同时允许用户选择任意合法非空组合。
        function configAllowsSelectableClassicMethodCombinations(testCase)
            cfg = config('smoke');
            expected = ["oracle", "nsb", "frozen_mvdr", "scm_mvdr", ...
                "nlms", "rls", "pi", "dl_mvdr", "esp", "wco", ...
                "static_incm_mvdr"];

            testCase.verifyEqual(cfg.experiments.classicJitterMethodNames, expected);
            testCase.verifyGreaterThan( ...
                cfg.experiments.periodCounts.jitterSweep, 0);
            testCase.verifyTrue(isfield(cfg.beamforming, 'classic'));

            subset = cfg;
            subset.experiments.classicJitterMethodNames = ...
                ["scm_mvdr", "wco", "oracle"];
            validated = uav.config.validate(subset);
            testCase.verifyEqual(validated.experiments.classicJitterMethodNames, ...
                ["scm_mvdr", "wco", "oracle"]);

            duplicate = cfg;
            duplicate.experiments.classicJitterMethodNames = ...
                ["oracle", "oracle"];
            testCase.verifyError(@() uav.config.validate(duplicate), ...
                'uav:config:InvalidExperimentParameter');

            unknown = cfg;
            unknown.experiments.classicJitterMethodNames = "unknown";
            testCase.verifyError(@() uav.config.validate(unknown), ...
                'uav:config:InvalidExperimentParameter');
        end

        % NSB必须在标称期望方向保持单位响应，并在标称干扰方向形成零陷。
        function nullSteeringSatisfiesNominalLinearConstraints(testCase)
            cfg = config('smoke');

            [weight, diagnostics] = uav.beamforming.nullSteering( ...
                cfg.derived.nominalSteering, cfg);
            responses = cfg.derived.nominalSteering' * weight;

            testCase.verifyEqual(responses(1), 1, 'AbsTol', 1e-10);
            testCase.verifyEqual(responses(2:end), ...
                zeros(cfg.derived.J, 1), 'AbsTol', 1e-10);
            testCase.verifyLessThan(diagnostics.constraintResidual, 1e-10);
        end

        % NLMS和RLS应使用已知训练参考逐快拍收敛，并返回可跨周期传递的状态。
        function adaptiveReferenceMethodsConvergeOnDeterministicTraining(testCase)
            cfg = config('smoke');
            M = cfg.derived.M;
            snapshots = repmat(eye(M), 1, 40);
            targetWeight = ((1:M).' + 1i * (M:-1:1).') / M;
            desiredReference = targetWeight' * snapshots;
            obs = struct('X', snapshots, ...
                'pilot', desiredReference / sqrt(cfg.derived.desiredPower));

            [nlmsWeight, nlmsState] = uav.beamforming.adaptiveReference( ...
                "nlms", obs, [], cfg);
            [rlsWeight, rlsState] = uav.beamforming.adaptiveReference( ...
                "rls", obs, [], cfg);

            nlmsError = mean(abs(desiredReference - nlmsWeight' * snapshots).^2);
            rlsError = mean(abs(desiredReference - rlsWeight' * snapshots).^2);
            testCase.verifyLessThan(nlmsError, 1e-6);
            testCase.verifyLessThan(rlsError, 1e-6);
            testCase.verifyEqual(nlmsState.method, "nlms");
            testCase.verifyEqual(rlsState.method, "rls");
            testCase.verifySize(rlsState.inverseCorrelation, [M, M]);
        end

        % PI仅固定配置参考阵元，不引入期望方向无失真约束。
        function powerInversionFixesConfiguredReferenceElement(testCase)
            cfg = config('smoke');
            M = cfg.derived.M;
            covariance = diag(1:M) + 0.05 * ones(M);

            [weight, diagnostics] = uav.beamforming.powerInversion( ...
                covariance, cfg);

            reference = cfg.beamforming.classic.pi.referenceElement;
            testCase.verifyEqual(weight(reference), 1, 'AbsTol', 1e-12);
            testCase.verifyEqual(diagnostics.referenceElement, reference);
        end

        % ESP修正后的导向矢量必须位于配置的主特征子空间内。
        function eigenspaceProjectionUsesConfiguredSignalSubspace(testCase)
            cfg = config('smoke');
            cfg.beamforming.classic.esp.signalSubspaceDimension = 2;
            cfg = uav.config.validate(cfg);
            M = cfg.derived.M;
            [basis, ~] = qr(complex(hilb(M), flipud(hilb(M))), 0);
            covariance = basis * diag(M:-1:1) * basis';
            reference = cfg.derived.nominalSteeringReference;

            [weight, diagnostics] = uav.beamforming.eigenspaceProjection( ...
                covariance, reference, cfg);

            residual = (eye(M) - diagnostics.signalSubspace * ...
                diagnostics.signalSubspace') * diagnostics.projectedReference;
            testCase.verifyLessThan(norm(residual), 1e-10);
            testCase.verifyEqual(weight' * diagnostics.projectedReference, ...
                1, 'AbsTol', 1e-9);
        end

        % 经典WCO不使用额外gamma正则项，并满足球形不确定集最坏响应约束。
        function wcoUsesZeroRegularizationAndSatisfiesWorstCaseResponse(testCase)
            cfg = config('smoke');
            cfg.beamforming.epsilonS = 0.2;
            reference = cfg.derived.nominalSteeringReference;

            [weight, diagnostics] = uav.beamforming.solve( ...
                "wco", eye(cfg.derived.M), reference, cfg);
            margin = real(weight' * reference) - ...
                cfg.beamforming.epsilonS * norm(weight);

            testCase.verifyEqual(diagnostics.method, "wco");
            testCase.verifyEqual(diagnostics.gamma, 0, 'AbsTol', 0);
            testCase.verifyGreaterThanOrEqual(margin, ...
                1 - cfg.beamforming.solver.feasibilityTolerance);
        end

        % 经典运行器必须只执行所选组合、保持顺序并让所有方法共享同一快拍。
        function classicRunnerUsesSelectedOrderAndSharedObservation(testCase)
            cfg = compactClassicConfig(testCase);
            scenario = uav.scenario.generate(cfg, ...
                uav.util.makeStream(cfg.random.baseSeed, 610001));
            obs = uav.scenario.generatePeriod(cfg, scenario, 1, ...
                uav.util.makeStream(cfg.random.baseSeed, 610101));
            methods = ["oracle", "nsb", "nlms", "static_incm_mvdr"];

            [record, nextState] = uav.experiment.runClassicMethods( ...
                obs, [], cfg, methods);

            testCase.verifyEqual(record.methodNames, methods);
            testCase.verifySize(record.sinrLinear, [1, numel(methods)]);
            testCase.verifyTrue(all(isfinite(record.sinrLinear)));
            testCase.verifyEqual(record.snapshotSignatures, ...
                repmat(record.snapshotSignatures(1, :), numel(methods), 1), ...
                'AbsTol', 0);
            testCase.verifyEqual(record.diagnostics{4}.beamforming.method, "mvdr");
            testCase.verifyEqual(record.diagnostics{4}.reconstruction.method, "static");
            testCase.verifyEqual(nextState.nlms.method, "nlms");
            testCase.verifyFalse(any(contains( ...
                string(fieldnames(nextState)), "truth")));
        end

        % 独立实验必须按配置子集生成统计量和完整结果文件，不加入固定五实验入口。
        function standaloneClassicJitterExperimentPublishesSelectedMethods(testCase)
            cfg = compactClassicConfig(testCase);
            cfg.experiments.jitterSweepDeg = [0, 1];
            cfg.experiments.classicJitterMethodNames = ...
                ["scm_mvdr", "dl_mvdr", "static_incm_mvdr"];
            cfg = uav.config.validate(cfg);

            result = exp06_classic_method_jitter_sweep(cfg);

            testCase.verifyEqual(result.methodNames, ...
                cfg.experiments.classicJitterMethodNames);
            testCase.verifySize(result.raw.sinrLinear, [1, 1, 2, 3]);
            testCase.verifySize(result.summary.sinr.medianDb, [2, 3]);
            artifactPaths = string(struct2cell(result.artifacts));
            artifactPaths = artifactPaths(endsWith(artifactPaths, ...
                [".mat"; ".csv"; ".fig"; ".pdf"; ".png"; ".json"]));
            testCase.verifyEqual(numel(artifactPaths), 6);
            testCase.verifyTrue(all(isfile(artifactPaths)));
            testCase.verifyFalse(any(uav.experiment.entryPoints() == ...
                "exp06_classic_method_jitter_sweep"));
        end
    end
end

function cfg = compactClassicConfig(testCase)
%COMPACTCLASSICCONFIG 构造经典方法测试专用的小规模配置和临时输出目录。
cfg = config('smoke');
cfg.timeline.K = 2;
cfg.timeline.L = 8;
cfg.monteCarlo.nTrials = 1;
cfg.spectrum.gridSize = 9;
cfg.beamforming.calibration.sampleCount = 16;
cfg.experiments.periodCounts.modelValidation = 2;
cfg.experiments.periodCounts.jitterSweep = 1;
cfg.experiments.periodCounts.snapshotSweep = 1;
cfg.experiments.periodCounts.dynamicStability = 2;
cfg.experiments.periodCounts.ablation = 1;
cfg.experiments.bootstrapReplicates = 12;
cfg.plot.visible = false;
cfg.plot.closeAfterExport = true;
cfg = uav.config.validate(cfg);
outputRoot = tempname;
mkdir(outputRoot);
testCase.addTeardown(@() cleanupClassicDirectory(outputRoot));
cfg.output.rootDir = outputRoot;
cfg.output.rawDir = fullfile(outputRoot, 'raw');
cfg.output.tableDir = fullfile(outputRoot, 'tables');
cfg.output.figureDir = fullfile(outputRoot, 'figures');
cfg.output.logDir = fullfile(outputRoot, 'logs');
cfg.output.cacheDir = fullfile(outputRoot, 'cache');
end

function cleanupClassicDirectory(pathValue)
%CLEANUPCLASSICDIRECTORY 删除测试独占且可完全重建的临时结果目录。
if isfolder(pathValue)
    rmdir(pathValue, 's');
end
end
