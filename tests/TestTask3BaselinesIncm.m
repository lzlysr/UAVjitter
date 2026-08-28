classdef TestTask3BaselinesIncm < matlab.unittest.TestCase
    %TESTTASK3BASELINESINCM 第三阶段：基线权值、Capon谱和静态INCM测试。
    % 既验证闭式公式，也覆盖病态SCM、空干扰、L<M和矩形UPA等边界。
    methods (TestMethodSetup)
        % 每个测试使用隔离路径，避免依赖用户MATLAB当前目录。
        function addSourcePath(testCase)
            rootDir = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(rootDir));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(rootDir, 'src')));
        end
    end

    methods (Test)
        % 配置必须提供静态重构所需的加载、噪声和结构修复参数。
        function configSuppliesStaticIncmNumericalParameters(testCase)
            cfg = config('smoke');

            testCase.verifyEqual(cfg.spectrum.noiseEigenvalueFraction, 0.5);
            testCase.verifyEqual(cfg.incm.noiseFloorRatio, 1e-6);
            testCase.verifyEqual(cfg.beamforming.minimumReciprocalCondition, 1e-8);
        end

        % 非法加载、比例和收缩参数必须在算法运行前被validate拒绝。
        function validationRejectsInvalidStaticIncmNumericalParameters(testCase)
            cfg = config('smoke');
            cfg.spectrum.noiseEigenvalueFraction = 0;
            testCase.verifyError(@() uav.config.validate(cfg), ...
                'uav:config:InvalidAlgorithmParameter');

            cfg = config('smoke');
            cfg.incm.noiseFloorRatio = -1;
            testCase.verifyError(@() uav.config.validate(cfg), ...
                'uav:config:InvalidAlgorithmParameter');

            cfg = config('smoke');
            cfg.beamforming.minimumReciprocalCondition = 1;
            testCase.verifyError(@() uav.config.validate(cfg), ...
                'uav:config:InvalidAlgorithmParameter');
        end

        % Capon网格必须为不小于3的奇数，确保包含中心节点且步长有定义。
        function validationRequiresOddCaponGridOfAtLeastThree(testCase)
            invalidGridSizes = [1, 2, 20];
            for gridSize = invalidGridSizes
                cfg = config('smoke');
                cfg.spectrum.gridSize = gridSize;
                testCase.verifyError(@() uav.config.validate(cfg), ...
                    'uav:config:InvalidPositiveInteger');
            end
        end

        % 标准MVDR结果应与手工闭式解一致，并满足无失真响应w^H*a=1。
        function mvdrHasClosedFormDistortionlessResponse(testCase)
            cfg = config('smoke');
            R = diag([2, 3]);
            aRef = [1; 1];

            [w, diagnostics] = uav.beamforming.solve("mvdr", R, aRef, cfg);

            testCase.verifyEqual(w, [3; 2] / 5, 'AbsTol', 1e-12);
            testCase.verifyEqual(w' * aRef, 1, 'AbsTol', 1e-12);
            testCase.verifyEqual(diagnostics.loading, 0, 'AbsTol', 0);
        end

        % 奇异SCM经DL-MVDR后应得到有限且满足响应约束的权值。
        function dlmvdrProducesFiniteDistortionlessWeightsForSingularScm(testCase)
            cfg = config('smoke');
            R = [1, 1; 1, 1];
            aRef = [1; -1];

            [w, diagnostics] = uav.beamforming.solve("dlmvdr", R, aRef, cfg);

            testCase.verifyTrue(all(isfinite(w)));
            testCase.verifyEqual(w' * aRef, 1, 'AbsTol', 1e-12);
            testCase.verifyGreaterThan(diagnostics.loading, 0);
        end

        % L<M造成秩亏时，标准MVDR仅施加必要数值加载且不发出求解警告。
        function mvdrNumericallyLoadsUnderdeterminedScmWithoutWarning(testCase)
            cfg = config('smoke');
            cfg.timeline.L = 3;
            cfg = uav.config.validate(cfg);
            scenario = uav.scenario.generate(cfg, uav.util.makeStream(53, 1));
            obs = uav.scenario.generatePeriod(cfg, scenario, 1, uav.util.makeStream(53, 2));
            aRef = cfg.derived.nominalSteeringReference;

            testCase.verifyWarningFree(@() uav.beamforming.solve("mvdr", obs.Rx, aRef, cfg));
            [w, diagnostics] = uav.beamforming.solve("mvdr", obs.Rx, aRef, cfg);

            testCase.verifyLessThan(cfg.timeline.L, cfg.derived.M);
            testCase.verifyGreaterThan(diagnostics.loading, 0);
            testCase.verifyGreaterThanOrEqual(diagnostics.reciprocalCondition, ...
                cfg.beamforming.minimumReciprocalCondition * (1 - 1e-8));
            testCase.verifyTrue(all(isfinite(w)));
            testCase.verifyLessThan(norm(w), 1 / sqrt(eps));
            testCase.verifyEqual(w' * aRef, 1, 'AbsTol', 1e-12);
        end

        % 用手算功率验证SINR线性代数，并确认J=0时干扰泄漏定义为0。
        function metricsUsesLinearSinrAlgebraAndHandlesZeroInterferers(testCase)
            cfg = config('smoke');
            w = [0.5; 2];
            truth = struct( ...
                'trueASVs', eye(2), ...
                'desiredPower', 4, ...
                'interfererPowers', 9, ...
                'trueIncm', diag([1, 10]));

            metrics = uav.metrics.evaluate(w, truth, cfg);

            testCase.verifyEqual(metrics.desiredGain, 0.25, 'AbsTol', 1e-12);
            testCase.verifyEqual(metrics.interferenceLeakage, 36, 'AbsTol', 1e-12);
            testCase.verifyEqual(metrics.sinrLinear, 1 / 40.25, 'AbsTol', 1e-12);
            testCase.verifyEqual(metrics.sinrDb, 10 * log10(1 / 40.25), 'AbsTol', 1e-12);

            noInterfererTruth = struct( ...
                'trueASVs', [1; 0], ...
                'desiredPower', 4, ...
                'interfererPowers', zeros(1, 0), ...
                'trueIncm', eye(2));
            noInterfererMetrics = uav.metrics.evaluate(w, noInterfererTruth, cfg);
            testCase.verifyEqual(noInterfererMetrics.interferenceLeakage, 0, 'AbsTol', 0);
            testCase.verifyEqual(noInterfererMetrics.sinrLinear, 1 / 4.25, 'AbsTol', 1e-12);
        end

        % 相同真值下Oracle MVDR输出SINR不应低于给定非Oracle权值。
        function oracleMvdrSinrIsNotBelowNominalNonOracleWeight(testCase)
            cfg = config('smoke');
            truth = struct( ...
                'trueASVs', [1, 0; 1, 1], ...
                'desiredPower', 10, ...
                'interfererPowers', 19, ...
                'trueIncm', diag([1, 20]));

            oracleWeight = uav.beamforming.solve("mvdr", truth.trueIncm, ...
                truth.trueASVs(:, 1), cfg);
            nominalWeight = uav.beamforming.solve("dlmvdr", diag([1, 20]), ...
                [1; 0], cfg);
            oracleMetrics = uav.metrics.evaluate(oracleWeight, truth, cfg);
            nominalMetrics = uav.metrics.evaluate(nominalWeight, truth, cfg);

            testCase.verifyGreaterThanOrEqual(oracleMetrics.sinrLinear, ...
                nominalMetrics.sinrLinear - 1e-12);
        end

        % 第三阶段非鲁棒算法不得受尚未使用的epsilonS字段影响。
        function taskThreeNonRobustAlgorithmsIgnoreEpsilonS(testCase)
            cfg = config('smoke');
            changedCfg = cfg;
            changedCfg.beamforming.epsilonS = 99;
            R = diag([2, 3]);
            aRef = [1; 1];
            scenario = uav.scenario.generate(cfg, uav.util.makeStream(67, 1));
            obs = uav.scenario.generatePeriod(cfg, scenario, 1, uav.util.makeStream(67, 2));

            mvdr = uav.beamforming.solve("mvdr", R, aRef, cfg);
            changedMvdr = uav.beamforming.solve("mvdr", R, aRef, changedCfg);
            dlmvdr = uav.beamforming.solve("dlmvdr", R, aRef, cfg);
            changedDlmvdr = uav.beamforming.solve("dlmvdr", R, aRef, changedCfg);
            staticIncm = uav.incm.reconstruct("static", obs, [], cfg);
            changedStaticIncm = uav.incm.reconstruct("static", obs, [], changedCfg);

            testCase.verifyEqual(changedMvdr, mvdr, 'AbsTol', 0);
            testCase.verifyEqual(changedDlmvdr, dlmvdr, 'AbsTol', 0);
            testCase.verifyEqual(changedStaticIncm, staticIncm, 'AbsTol', 0);
        end

        % 合成单源协方差的Capon峰必须落在已知二维方向余弦附近。
        function loadedCaponLocalizesKnownSourceOnVisibleGrid(testCase)
            cfg = config('smoke');
            cfg.array.Nx = 2;
            cfg.array.Ny = 2;
            cfg.spectrum.gridSize = 21;
            cfg = uav.config.validate(cfg);
            targetXi = [0.2; -0.5];
            a = uav.array.steeringUPA(targetXi, cfg.array.Nx, cfg.array.Ny, ...
                cfg.array.spacingOverLambda);
            Rx = 100 * (a * a') + eye(cfg.derived.M);

            [powerMap, diagnostics] = uav.spectrum.loadedCapon(Rx, cfg);

            [~, linearIndex] = max(powerMap(:), [], 'omitnan');
            [row, column] = ind2sub(size(powerMap), linearIndex);
            peakXi = [diagnostics.xiX(row, column); diagnostics.xiY(row, column)];
            testCase.verifyEqual(peakXi, targetXi, 'AbsTol', 1e-12);
            testCase.verifyTrue(diagnostics.visibleMask(row, column));
        end

        % 矩形UPA二维谱在可见圆盘内有限，圆盘外必须统一为NaN。
        function loadedCaponMasksOutsideDiskForRectangularUpa(testCase)
            cfg = config('smoke');
            cfg.array.Nx = 2;
            cfg.array.Ny = 3;
            cfg.spectrum.gridSize = 11;
            cfg = uav.config.validate(cfg);

            [powerMap, diagnostics] = uav.spectrum.loadedCapon(eye(cfg.derived.M), cfg);

            testCase.verifySize(powerMap, [11, 11]);
            testCase.verifyTrue(any(~diagnostics.visibleMask, 'all'));
            testCase.verifyTrue(all(isnan(powerMap(~diagnostics.visibleMask))));
            testCase.verifyTrue(all(isfinite(powerMap(diagnostics.visibleMask))));
        end

        % 缓存流形和批量Capon必须与逐点参考实现数值一致。
        function cachedGridManifoldAndBatchCaponMatchPointwiseReference(testCase)
            cfg = config('smoke');
            cfg.array.Nx = 2;
            cfg.array.Ny = 3;
            cfg.spectrum.gridSize = 11;
            cfg = uav.config.validate(cfg);
            grid = cfg.derived.spectrumGrid;

            testCase.verifySize(grid.manifold, ...
                [cfg.derived.M, nnz(grid.visibleMask)]);
            testCase.verifyEqual(grid.directionCosines, ...
                [grid.xiX(grid.visibleMask).'; grid.xiY(grid.visibleMask).'], ...
                'AbsTol', 0);
            pointwiseManifold = zeros(size(grid.manifold));
            for index = 1:size(grid.directionCosines, 2)
                pointwiseManifold(:, index) = uav.array.steeringUPA( ...
                    grid.directionCosines(:, index), cfg.array.Nx, ...
                    cfg.array.Ny, cfg.array.spacingOverLambda);
            end
            testCase.verifyEqual(grid.manifold, pointwiseManifold, 'AbsTol', 0);

            seedMatrix = reshape(1:cfg.derived.M^2, cfg.derived.M, []);
            Rx = seedMatrix * seedMatrix' + eye(cfg.derived.M);
            [powerMap, diagnostics] = uav.spectrum.loadedCapon(Rx, cfg);
            reference = nan(size(powerMap));
            for index = 1:numel(grid.visibleLinearIndices)
                a = pointwiseManifold(:, index);
                reference(grid.visibleLinearIndices(index)) = ...
                    1 / real(a' * (diagnostics.loadedCovariance \ a));
            end
            testCase.verifyEqual(powerMap, reference, 'RelTol', 1e-12, ...
                'AbsTol', 1e-12);
            testCase.verifyTrue(diagnostics.usedCachedManifold);
        end

        % 多个J/L/Nx/Ny边界下静态INCM必须有限、Hermitian且半正定。
        function staticIncmIsHermitianPsdAndFiniteForBoundaryCases(testCase)
            cases = { ...
                struct('Nx', 3, 'Ny', 3, 'L', 3, 'directions', [40, -25], 'inr', 20); ...
                struct('Nx', 3, 'Ny', 3, 'L', 3, 'directions', zeros(0, 2), 'inr', []); ...
                struct('Nx', 1, 'Ny', 2, 'L', 3, ...
                    'directions', [40, -25; -55, -20; 70, -30], 'inr', [20, 15, 12])};

            for caseIndex = 1:numel(cases)
                item = cases{caseIndex};
                cfg = config('smoke');
                cfg.array.Nx = item.Nx;
                cfg.array.Ny = item.Ny;
                cfg.timeline.L = item.L;
                cfg.sources.interfererAzElDeg = item.directions;
                cfg.sources.interfererInrDb = item.inr;
                cfg = uav.config.validate(cfg);
                scenario = uav.scenario.generate(cfg, uav.util.makeStream(71, caseIndex));
                obs = uav.scenario.generatePeriod(cfg, scenario, 1, ...
                    uav.util.makeStream(71, caseIndex + 10));

                [RinHat, state, diagnostics] = uav.incm.reconstruct("static", obs, [], cfg);

                testCase.verifyEqual(RinHat, RinHat', 'AbsTol', 1e-12);
                testCase.verifyGreaterThanOrEqual(min(eig(RinHat)), ...
                    diagnostics.eigenvalueFloor * (1 - 1e-10));
                testCase.verifyTrue(all(isfinite(RinHat), 'all'));
                testCase.verifyEqual(state.method, "static");
            end
        end

        % L<M时零特征值不含噪声信息，应使用配置接收机噪声先验。
        function underdeterminedNoiseUsesConfiguredPrior(testCase)
            cfg = config('smoke');
            cfg.timeline.L = 3;
            cfg.sources.interfererAzElDeg = zeros(0, 2);
            cfg.sources.interfererInrDb = [];
            cfg = uav.config.validate(cfg);
            X = sqrt(cfg.derived.M * cfg.derived.noisePower) * ...
                [eye(cfg.timeline.L); zeros(cfg.derived.M - cfg.timeline.L, cfg.timeline.L)];
            obs = struct('Rx', X * X' / cfg.timeline.L, 'X', X);

            [~, ~, diagnostics] = uav.incm.reconstruct("static", obs, [], cfg);

            testCase.verifyEqual(diagnostics.noiseEstimate, cfg.derived.noisePower, ...
                'AbsTol', 1e-12);
            testCase.verifyGreaterThan(diagnostics.noiseEstimate, diagnostics.noiseFloor);
            testCase.verifyEqual(diagnostics.effectiveSnapshotCount, cfg.timeline.L);
            testCase.verifyEqual(diagnostics.noiseEstimateSource, ...
                "configured-underdetermined");
        end

        % 强干扰且L<M时仍不能把非零小特征值误当成可靠噪声估计。
        function underdeterminedStrongInterferenceUsesConfiguredNoisePrior(testCase)
            for snapshotCount = [2, 4, 8]
                cfg = config('smoke');
                cfg.timeline.L = snapshotCount;
                cfg = uav.config.validate(cfg);
                scenario = uav.scenario.generate(cfg, uav.util.makeStream(73, snapshotCount));
                obs = uav.scenario.generatePeriod(cfg, scenario, 1, ...
                    uav.util.makeStream(74, snapshotCount));

                [~, ~, diagnostics] = uav.incm.reconstruct("static", obs, [], cfg);

                testCase.verifyEqual(diagnostics.noiseEstimate, cfg.derived.noisePower, ...
                    'AbsTol', 0);
                testCase.verifyEqual(diagnostics.noiseEstimateSource, ...
                    "configured-underdetermined");
            end
        end

        % L>=M的满秩纯噪声SCM应启用较小特征值估计器。
        function fullRankNoiseScmUsesLowerEigenvalueEstimator(testCase)
            cfg = config('smoke');
            cfg.timeline.L = 12;
            cfg.sources.interfererAzElDeg = zeros(0, 2);
            cfg.sources.interfererInrDb = [];
            cfg = uav.config.validate(cfg);
            X = sqrt(cfg.timeline.L * cfg.derived.noisePower) * ...
                [eye(cfg.derived.M), zeros(cfg.derived.M, cfg.timeline.L - cfg.derived.M)];
            obs = struct('Rx', X * X' / cfg.timeline.L, 'X', X);

            [~, ~, diagnostics] = uav.incm.reconstruct("static", obs, [], cfg);

            testCase.verifyEqual(diagnostics.noiseEstimate, cfg.derived.noisePower, ...
                'AbsTol', 1e-12);
            testCase.verifyEqual(diagnostics.noiseEstimateSource, "lower-eigenvalues");
        end

        % 原始快拍阵元行数与M不一致时必须拒绝，避免Rx与X语义分裂。
        function staticIncmRejectsObservationWithMismatchedSnapshotRows(testCase)
            cfg = config('smoke');
            obs = struct('Rx', eye(cfg.derived.M), ...
                'X', zeros(cfg.derived.M - 1, cfg.timeline.L));

            testCase.verifyError(@() uav.incm.reconstruct("static", obs, [], cfg), ...
                'uav:incm:InvalidObservation');
        end

        % 精确白噪声协方差下，扣除Capon噪声底座后残差谱应为零。
        function staticIncmSubtractsCaponNoisePedestalForExactNoiseCovariance(testCase)
            cfg = config('smoke');
            cfg.sources.interfererAzElDeg = zeros(0, 2);
            cfg.sources.interfererInrDb = [];
            cfg = uav.config.validate(cfg);
            noiseOnlyObservation = struct('Rx', cfg.derived.noisePower * eye(cfg.derived.M));

            [RinHat, ~, diagnostics] = uav.incm.reconstruct("static", ...
                noiseOnlyObservation, [], cfg);

            testCase.verifyEqual(RinHat, cfg.derived.noisePower * eye(cfg.derived.M), ...
                'AbsTol', 1e-12);
            testCase.verifyEqual(diagnostics.residualPowerMap(diagnostics.spectrum.visibleMask), ...
                zeros(nnz(diagnostics.spectrum.visibleMask), 1), 'AbsTol', 1e-12);
            testCase.verifyEqual(RinHat - diag(diag(RinHat)), zeros(cfg.derived.M), ...
                'AbsTol', 1e-12);
            testCase.verifyEqual(diagnostics.noiseEstimateSource, ...
                "configured-missing-snapshots");
        end

        % 删除truth或把truth替换为毒化值都不应改变静态重构结果。
        function staticIncmDoesNotReadRemovedOrPoisonedTruth(testCase)
            cfg = config('smoke');
            scenario = uav.scenario.generate(cfg, uav.util.makeStream(79, 1));
            obs = uav.scenario.generatePeriod(cfg, scenario, 1, uav.util.makeStream(79, 2));
            truthRemoved = rmfield(obs, 'truth');
            truthPoisoned = obs;
            truthPoisoned.truth = struct('trueASVs', nan(cfg.derived.M, cfg.derived.J + 1), ...
                'trueIncm', nan(cfg.derived.M), 'desiredPower', -1, ...
                'interfererPowers', inf(1, cfg.derived.J));

            reference = uav.incm.reconstruct("static", obs, [], cfg);
            removed = uav.incm.reconstruct("static", truthRemoved, [], cfg);
            poisoned = uav.incm.reconstruct("static", truthPoisoned, [], cfg);

            testCase.verifyEqual(removed, reference, 'AbsTol', 0);
            testCase.verifyEqual(poisoned, reference, 'AbsTol', 0);
        end

        % 期望方向排除区不得进入均匀求积，面积权重应与网格单元一致。
        function staticIncmExcludesNominalSoiWithUniformGridQuadrature(testCase)
            cfg = config('smoke');
            scenario = uav.scenario.generate(cfg, uav.util.makeStream(81, 1));
            obs = uav.scenario.generatePeriod(cfg, scenario, 1, uav.util.makeStream(81, 2));

            [~, ~, diagnostics] = uav.incm.reconstruct("static", obs, [], cfg);

            desiredXi = cfg.derived.nominalDirectionCosines(:, 1);
            distanceSquared = (diagnostics.spectrum.xiX - desiredXi(1)).^2 + ...
                (diagnostics.spectrum.xiY - desiredXi(2)).^2;
            [~, nearestIndex] = min(distanceSquared(:));
            expectedCellArea = (2 / (cfg.spectrum.gridSize - 1))^2;
            testCase.verifyTrue(diagnostics.spectrum.visibleMask(nearestIndex));
            testCase.verifyTrue(diagnostics.excludedSoiMask(nearestIndex));
            testCase.verifyFalse(diagnostics.includedMask(nearestIndex));
            testCase.verifyEqual(diagnostics.spectrum.cellArea, expectedCellArea, ...
                'AbsTol', 1e-12);
        end

        % 相同输入静态重构逐元素复现，并只在提供估计INCM时报告NMSE。
        function staticIncmIsDeterministicAndReportsOptionalNmse(testCase)
            cfg = config('smoke');
            scenario = uav.scenario.generate(cfg, uav.util.makeStream(83, 1));
            obs = uav.scenario.generatePeriod(cfg, scenario, 1, uav.util.makeStream(83, 2));

            [first, ~, firstDiagnostics] = uav.incm.reconstruct("static", obs, [], cfg);
            [second, ~, secondDiagnostics] = uav.incm.reconstruct("static", obs, [], cfg);
            metricsTruth = obs.truth;
            metricsTruth.estimatedIncm = first;
            metrics = uav.metrics.evaluate(cfg.derived.nominalSteeringReference, metricsTruth, cfg);

            testCase.verifyEqual(first, second, 'AbsTol', 0);
            testCase.verifyEqual(firstDiagnostics.noiseEstimate, secondDiagnostics.noiseEstimate, ...
                'AbsTol', 0);
            testCase.verifyTrue(isfield(metrics, 'incmNmse'));
            testCase.verifyGreaterThanOrEqual(metrics.incmNmse, 0);
        end
    end
end
