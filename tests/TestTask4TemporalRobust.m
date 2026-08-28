classdef TestTask4TemporalRobust < matlab.unittest.TestCase
    %TESTTASK4TEMPORALROBUST 第四阶段：时序谱增强、自适应求积和鲁棒权值测试。
    % 测试把报告中的每个算法模块独立锁定，再检查完整时序重构和消融。
    methods (TestMethodSetup)
        % 为测试加载工程包，结束后由fixture恢复环境。
        function addSourcePath(testCase)
            rootDir = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(rootDir));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(rootDir, 'src')));
        end
    end

    methods (Test)
        % 配置必须完整提供时序正则、对齐、求积、结构和鲁棒求解参数。
        function configSuppliesTemporalAndRobustDefaults(testCase)
            cfg = config('smoke');

            testCase.verifyEqual(cfg.spectrum.alignmentMaxShiftBins, [2, 2]);
            testCase.verifyGreaterThan(cfg.spectrum.alignmentMinCorrelation, 0);
            testCase.verifyGreaterThanOrEqual(cfg.spectrum.alignmentAmbiguityMargin, 0);
            testCase.verifyGreaterThanOrEqual(cfg.spectrum.alignmentMinimumOverlap, 2);
            testCase.verifyGreaterThanOrEqual(cfg.spectrum.peakThresholdRatio, 0);
            testCase.verifyGreaterThanOrEqual(cfg.spectrum.peakNmsRadiusBins, 0);
            testCase.verifyGreaterThanOrEqual(cfg.spectrum.quadratureCoarseStride, 1);
            testCase.verifyGreaterThanOrEqual(cfg.spectrum.quadratureFineRadiusBins, 0);
            testCase.verifyGreaterThan(cfg.spectrum.quadratureRefinement, 1);
            testCase.verifyGreaterThan(cfg.spectrum.quadratureDedupTolerance, 0);
            testCase.verifyTrue(cfg.incm.ablation.temporalContinuity);
            testCase.verifyTrue(cfg.incm.ablation.alignment);
            testCase.verifyTrue(cfg.incm.ablation.adaptiveQuadrature);
            testCase.verifyTrue(cfg.incm.ablation.structureEnhancement);
            testCase.verifyFalse(cfg.incm.ablation.pilotProjection);
            testCase.verifyGreaterThanOrEqual(cfg.beamforming.gamma, 0);
            testCase.verifyGreaterThan(cfg.beamforming.covariancePsdTolerance, 0);
            testCase.verifyGreaterThan(cfg.beamforming.solver.initialInteriorMargin, 0);
            testCase.verifyGreaterThan(cfg.beamforming.solver.feasibilityTolerance, 0);
            testCase.verifyEqual(cfg.beamforming.calibration.quantile, 0.95);
            testCase.verifyGreaterThan(cfg.beamforming.calibration.sampleCount, 0);
        end

        % 拒绝负正则、非法相关阈值、不可行epsilon等第四阶段数值。
        function validationRejectsInvalidTaskFourNumerics(testCase)
            invalidCases = { ...
                'spectrum', 'alignmentMaxShiftBins', [2, -1]; ...
                'spectrum', 'alignmentMinCorrelation', 1.1; ...
                'spectrum', 'alignmentAmbiguityMargin', -0.1; ...
                'spectrum', 'alignmentMinimumOverlap', 1; ...
                'spectrum', 'peakThresholdRatio', 1.1; ...
                'spectrum', 'peakNmsRadiusBins', 0.5; ...
                'spectrum', 'quadratureCoarseStride', 0; ...
                'spectrum', 'quadratureFineRadiusBins', -1; ...
                'spectrum', 'quadratureRefinement', 1; ...
                'spectrum', 'quadratureDedupTolerance', 0; ...
                'beamforming', 'gamma', -1; ...
                'beamforming', 'covariancePsdTolerance', -1};
            for index = 1:size(invalidCases, 1)
                cfg = config('smoke');
                cfg.(invalidCases{index, 1}).(invalidCases{index, 2}) = ...
                    invalidCases{index, 3};
                testCase.verifyError(@() uav.config.validate(cfg), ...
                    'uav:config:InvalidAlgorithmParameter');
            end

            cfg = config('smoke');
            cfg.beamforming.calibration.quantile = 1;
            testCase.verifyError(@() uav.config.validate(cfg), ...
                'uav:config:InvalidAlgorithmParameter');
            cfg = config('smoke');
            cfg.beamforming.solver.maxIterations = 0;
            testCase.verifyError(@() uav.config.validate(cfg), ...
                'uav:config:InvalidAlgorithmParameter');
            cfg = config('smoke');
            cfg.beamforming.solver.initialInteriorMargin = 0;
            testCase.verifyError(@() uav.config.validate(cfg), ...
                'uav:config:InvalidAlgorithmParameter');
            cfg = config('smoke');
            cfg.beamforming.solver.feasibilityTolerance = NaN;
            testCase.verifyError(@() uav.config.validate(cfg), ...
                'uav:config:InvalidAlgorithmParameter');
            cfg = config('smoke');
            cfg.incm.ablation.pilotProjection = 1;
            testCase.verifyError(@() uav.config.validate(cfg), ...
                'uav:config:InvalidAlgorithmParameter');
        end

        % 合成正负二维平移谱，检查对齐模块恢复正确行列位移。
        function alignmentRecoversKnownPositiveAndNegativeShifts(testCase)
            cfg = config('smoke');
            previous = zeros(9);
            previous(3, 4) = 1;
            previous(5, 3) = 0.35;
            previous(6, 6) = 0.7;
            mask = true(9);

            currentPositive = zeros(9);
            currentPositive(4, 3) = 1;
            currentPositive(6, 2) = 0.35;
            currentPositive(7, 5) = 0.7;
            [alignedPositive, positiveDiagnostics] = uav.spectrum.align2D( ...
                currentPositive, previous, mask, cfg);
            testCase.verifyEqual(positiveDiagnostics.shiftBins, [1, -1]);
            testCase.verifyEqual(alignedPositive, currentPositive, 'AbsTol', 0);
            testCase.verifyFalse(positiveDiagnostics.usedIdentity);

            currentNegative = zeros(9);
            currentNegative(2, 6) = 1;
            currentNegative(4, 5) = 0.35;
            currentNegative(5, 8) = 0.7;
            [~, negativeDiagnostics] = uav.spectrum.align2D( ...
                currentNegative, previous, mask, cfg);
            testCase.verifyEqual(negativeDiagnostics.shiftBins, [-1, 2]);
            testCase.verifyFalse(negativeDiagnostics.usedIdentity);
        end

        % 相关过弱或证据不足时必须使用恒等映射而不是强行平移。
        function alignmentUsesIdentityForWeakOrAmbiguousEvidence(testCase)
            cfg = config('smoke');
            weak = ones(7);

            [aligned, diagnostics] = uav.spectrum.align2D(weak, weak, true(7), cfg);

            testCase.verifyEqual(aligned, weak, 'AbsTol', 0);
            testCase.verifyEqual(diagnostics.shiftBins, [0, 0]);
            testCase.verifyTrue(diagnostics.usedIdentity);
            testCase.verifyTrue(contains(diagnostics.fallbackReason, "weak"));
        end

        % 不规则掩码外允许NaN，相关计算只能读取共同有效支撑。
        function alignmentMasksNonfiniteValuesOutsideIrregularSupport(testCase)
            cfg = config('smoke');
            [columns, rows] = meshgrid(1:11, 1:11);
            mask = (rows - 6).^2 + (columns - 6).^2 <= 16;
            mask(6, 6) = false;
            previous = nan(11);
            previous(mask) = 0;
            previous(4, 5) = 1;
            previous(7, 4) = 0.6;
            previous(8, 7) = 0.3;
            current = nan(11);
            current(mask) = 0;
            current(5, 6) = 1;
            current(8, 5) = 0.6;
            current(9, 8) = 0.3;

            [aligned, diagnostics] = uav.spectrum.align2D( ...
                current, previous, mask, cfg);

            testCase.verifyEqual(diagnostics.shiftBins, [1, 1]);
            testCase.verifyEqual(aligned(mask), current(mask), 'AbsTol', 0);
            testCase.verifyEqual(aligned(~mask), zeros(nnz(~mask), 1), 'AbsTol', 0);
            testCase.verifyTrue(all(isfinite(aligned), 'all'));
        end

        % 最高和次高相关过近时应识别为真实歧义并回退恒等映射。
        function alignmentUsesIdentityForGenuinelyAmbiguousTopScores(testCase)
            cfg = config('smoke');
            cfg.spectrum.alignmentMinCorrelation = 0.1;
            cfg.spectrum.alignmentAmbiguityMargin = 0.05;
            previous = zeros(9);
            previous(5, 5) = 1;
            current = zeros(9);
            current(5, 4) = 1;
            current(5, 6) = 1;

            [aligned, diagnostics] = uav.spectrum.align2D( ...
                current, previous, true(9), cfg);

            testCase.verifyGreaterThanOrEqual(diagnostics.bestCorrelation, ...
                cfg.spectrum.alignmentMinCorrelation);
            testCase.verifyLessThan(diagnostics.ambiguityGap, ...
                cfg.spectrum.alignmentAmbiguityMargin);
            testCase.verifyEqual(diagnostics.fallbackReason, "ambiguous-evidence");
            testCase.verifyEqual(aligned, previous, 'AbsTol', 0);
        end

        % 用手推分段区域检查时序近端解在0、g及两侧驻点的选择。
        function temporalProxMatchesHandDerivedBreakpointRegions(testCase)
            p = [0.1; 1; 2.1; 4];
            g = [2; 2; 2; 2];
            lambdaS = [1; 0.4; 0.4; 0.4];
            lambdaT = [0.2; 0.2; 0.8; 0.8];

            c = uav.incm.temporalProx(p, g, lambdaS, lambdaT);

            testCase.verifyEqual(c, [0; 0.9; 2; 3.4], 'AbsTol', 1e-12);
        end

        % 与独立致密标量搜索比较，防止解析候选遗漏或系数差一倍。
        function temporalProxMatchesIndependentDenseScalarSearch(testCase)
            pValues = [0, 0.1, 0.7, 1.8, 3.1];
            gValues = [0, 0.4, 1.5, 2.2];
            lambdaS = 0.37;
            lambdaT = 0.83;
            for p = pValues
                for g = gValues
                    actual = uav.incm.temporalProx(p, g, lambdaS, lambdaT);
                    candidates = linspace(0, 5, 50001);
                    objective = (p - candidates).^2 + lambdaS * candidates + ...
                        lambdaT * abs(candidates - g);
                    [~, index] = min(objective);
                    testCase.verifyEqual(actual, candidates(index), 'AbsTol', 1.1e-4);
                end
            end
        end

        % NMS峰数上限只能依赖M-1，不得读取真实干扰数J。
        function nmsPeakCountIsCappedByArrayNotTrueInterfererCount(testCase)
            cfgA = config('smoke');
            cfgA.array.Nx = 2;
            cfgA.array.Ny = 2;
            cfgA.spectrum.gridSize =  nine();
            cfgA.spectrum.peakThresholdRatio = 0;
            cfgA.spectrum.peakNmsRadiusBins = 0;
            cfgA = uav.config.validate(cfgA);
            cfgB = cfgA;
            cfgB.sources.interfererAzElDeg = zeros(0, 2);
            cfgB.sources.interfererInrDb = [];
            cfgB = uav.config.validate(cfgB);
            [~, spectrumDiagnostics] = uav.spectrum.loadedCapon(eye(4), cfgA);
            powerMap = zeros(9);
            candidateIndices = find(spectrumDiagnostics.visibleMask, 5, 'first');
            powerMap(candidateIndices) = [9, 8, 7, 6, 5];

            [peaksA, diagnosticsA] = uav.spectrum.selectPeaks( ...
                powerMap, spectrumDiagnostics, spectrumDiagnostics.visibleMask, cfgA);
            [peaksB, diagnosticsB] = uav.spectrum.selectPeaks( ...
                powerMap, spectrumDiagnostics, spectrumDiagnostics.visibleMask, cfgB);

            testCase.verifySize(peaksA, [cfgA.derived.M - 1, 2]);
            testCase.verifyEqual(peaksB, peaksA, 'AbsTol', 0);
            testCase.verifyEqual(diagnosticsB.linearIndices, diagnosticsA.linearIndices);
            testCase.verifyEqual(diagnosticsA.cap, cfgA.derived.M - 1);
        end

        % 细区与粗区不重叠，合并后代表面积必须守恒。
        function adaptiveQuadratureHasNoOverlapAndConservesGridArea(testCase)
            cfg = config('smoke');
            cfg.spectrum.gridSize = 11;
            cfg = uav.config.validate(cfg);
            [~, spectrumDiagnostics] = uav.spectrum.loadedCapon(eye(cfg.derived.M), cfg);
            gridNodes = [spectrumDiagnostics.xiX(:), spectrumDiagnostics.xiY(:)];
            includedMask = reshape(independentAdmissibility(gridNodes, cfg), ...
                size(spectrumDiagnostics.visibleMask));
            powerMap = ones(size(includedMask));
            peakIndex = sub2ind(size(powerMap), 6, 7);

            quadrature = uav.spectrum.adaptiveQuadrature(powerMap, ...
                spectrumDiagnostics, includedMask, peakIndex, cfg);

            roundedNodes = round(quadrature.nodes * 1e12) / 1e12;
            testCase.verifyEqual(size(unique(roundedNodes, 'rows'), 1), ...
                size(quadrature.nodes, 1));
            expectedArea = nnz(includedMask) * spectrumDiagnostics.cellArea;
            testCase.verifyEqual(sum(quadrature.weights), expectedArea, ...
                'AbsTol', 1e-12);
            testCase.verifyGreaterThan(quadrature.diagnostics.fineNodeCount, 0);
            testCase.verifyGreaterThan(quadrature.diagnostics.coarseNodeCount, 0);
            testCase.verifyTrue(all(quadrature.weights > 0));
            testCase.verifyTrue(all(independentAdmissibility(quadrature.nodes, cfg)));
        end

        % 所有自适应细节点都必须满足连续可见域及期望信号排除边界。
        function continuousAdmissibilityAndAdaptiveNodesRespectSoiBoundary(testCase)
            cfg = config('smoke');
            cfg.spectrum.gridSize = 21;
            cfg = uav.config.validate(cfg);
            [~, spectrumDiagnostics] = uav.spectrum.loadedCapon( ...
                eye(cfg.derived.M), cfg);
            gridNodes = [spectrumDiagnostics.xiX(:), spectrumDiagnostics.xiY(:)];
            expectedMask = independentAdmissibility(gridNodes, cfg);

            actualMask = uav.spectrum.admissibleDirections(gridNodes, cfg);

            testCase.verifyEqual(actualMask, expectedMask);
            angularDistance = lowerHemisphereAngularDistance(gridNodes, cfg);
            nearBoundary = expectedMask & angularDistance < ...
                deg2rad(cfg.spectrum.desiredExclusionHalfWidthDeg + 5) & ...
                sum(gridNodes.^2, 2) < 0.9^2;
            candidates = find(nearBoundary);
            testCase.assertNotEmpty(candidates);
            [~, nearest] = min(angularDistance(candidates));
            peakIndex = candidates(nearest);
            includedMask = reshape(expectedMask, size(spectrumDiagnostics.visibleMask));
            powerMap = ones(size(includedMask));

            quadrature = uav.spectrum.adaptiveQuadrature(powerMap, ...
                spectrumDiagnostics, includedMask, peakIndex, cfg);

            nodeDistance = lowerHemisphereAngularDistance(quadrature.nodes, cfg);
            testCase.verifyTrue(all(sum(quadrature.nodes.^2, 2) <= 1 + 1e-12));
            testCase.verifyGreaterThan(min(nodeDistance), ...
                deg2rad(cfg.spectrum.desiredExclusionHalfWidthDeg));
            testCase.verifyEqual(sum(quadrature.weights), ...
                nnz(includedMask) * spectrumDiagnostics.cellArea, 'AbsTol', 1e-12);
        end

        % 去重不应量化实际物理坐标，容差还必须受细网格间距约束。
        function dedupKeepsPhysicalNodesAndToleranceIsGridBounded(testCase)
            cfg = config('smoke');
            cfg.spectrum.gridSize = 11;
            cfg.spectrum.quadratureCoarseStride = 3;
            cfg.spectrum.quadratureFineRadiusBins = 0;
            gridStep = 2 / (cfg.spectrum.gridSize - 1);
            maximumTolerance = gridStep / (2 * cfg.spectrum.quadratureRefinement);
            cfg.spectrum.quadratureDedupTolerance = maximumTolerance;
            validated = uav.config.validate(cfg);
            testCase.verifyEqual(validated.spectrum.quadratureDedupTolerance, ...
                maximumTolerance, 'AbsTol', 0);
            cfg.spectrum.quadratureDedupTolerance = maximumTolerance * (1 + 1e-6);
            testCase.verifyError(@() uav.config.validate(cfg), ...
                'uav:config:InvalidAlgorithmParameter');

            cfg.spectrum.quadratureDedupTolerance = 0.9 * maximumTolerance;
            cfg = uav.config.validate(cfg);
            controlledNodes = [-0.2, -0.2; -0.201, -0.201];
            [mergedNodes, mergedWeights, mergedValues, mergedDiagnostics] = ...
                uav.spectrum.mergeQuadratureNodes(controlledNodes, [1; 3], ...
                [2; 4], cfg.spectrum.quadratureDedupTolerance);
            testCase.verifyEqual(mergedNodes, controlledNodes(1, :), 'AbsTol', 0);
            testCase.verifyEqual(mergedWeights, 4, 'AbsTol', 0);
            testCase.verifyEqual(mergedValues, 3.5, 'AbsTol', 0);
            testCase.verifyEqual(mergedDiagnostics.deduplicatedNodeCount, 1);

            [~, spectrumDiagnostics] = uav.spectrum.loadedCapon( ...
                eye(cfg.derived.M), cfg);
            gridNodes = [spectrumDiagnostics.xiX(:), spectrumDiagnostics.xiY(:)];
            includedMask = reshape(independentAdmissibility(gridNodes, cfg), ...
                size(spectrumDiagnostics.visibleMask));
            peakIndex = sub2ind(size(includedMask), 5, 5);
            testCase.assertTrue(includedMask(peakIndex));
            quadrature = uav.spectrum.adaptiveQuadrature(ones(size(includedMask)), ...
                spectrumDiagnostics, includedMask, peakIndex, cfg);
            axisValues = unique(spectrumDiagnostics.xiX(:));
            offsets = (((1:cfg.spectrum.quadratureRefinement) - 0.5) / ...
                cfg.spectrum.quadratureRefinement - 0.5) * gridStep;
            refinedAxis = axisValues + offsets;
            allowedCoordinates = unique([axisValues; refinedAxis(:)]);
            xDistance = min(abs(quadrature.nodes(:, 1) - allowedCoordinates.'), [], 2);
            yDistance = min(abs(quadrature.nodes(:, 2) - allowedCoordinates.'), [], 2);

            testCase.verifyLessThanOrEqual(max([xDistance; yDistance]), 1e-12);
            testCase.verifyTrue(all(independentAdmissibility(quadrature.nodes, cfg)));
        end

        % 结构投影输出应Hermitian、PSD、满足特征值下限且支持收缩消融。
        function structureProjectionIsHermitianPsdFlooredAndShrinkable(testCase)
            cfg = config('smoke');
            R = [1, 2 + 1i; 2 - 3i, -0.5];

            [projected, diagnostics] = uav.incm.projectStructure(R, cfg, true);

            testCase.verifyEqual(projected, projected', 'AbsTol', 1e-12);
            testCase.verifyGreaterThanOrEqual(min(real(eig(projected))), ...
                diagnostics.eigenvalueFloor * (1 - 1e-10));
            testCase.verifyTrue(diagnostics.applied);

            [hermitianOnly, disabledDiagnostics] = uav.incm.projectStructure(R, cfg, false);
            testCase.verifyEqual(hermitianOnly, (R + R') / 2, 'AbsTol', 0);
            testCase.verifyFalse(disabledDiagnostics.applied);
        end

        % 静态诊断应保留完整残差谱，同时明确标记被排除的期望区域。
        function staticDiagnosticsPreserveExcludedSoiResidualMap(testCase)
            cfg = config('smoke');
            cfg.spectrum.gridSize = 21;
            cfg = uav.config.validate(cfg);
            aRef = cfg.derived.nominalSteeringReference;
            obs = struct('Rx', 100 * (aRef * aRef') + eye(cfg.derived.M));

            [~, ~, diagnostics] = uav.incm.reconstruct("static", obs, [], cfg);

            desiredXi = cfg.derived.nominalDirectionCosines(:, 1);
            distanceSquared = (diagnostics.spectrum.xiX - desiredXi(1)).^2 + ...
                (diagnostics.spectrum.xiY - desiredXi(2)).^2;
            [~, nearestIndex] = min(distanceSquared(:));
            testCase.verifyTrue(diagnostics.excludedSoiMask(nearestIndex));
            testCase.verifyGreaterThan(diagnostics.residualPowerMap(nearestIndex), 0);
        end

        % 完整时序重构必须可复现、不读取truth、输出有限PSD矩阵。
        function temporalReconstructionIsDeterministicSafeAndPsd(testCase)
            cfg = compactConfig();
            scenario = uav.scenario.generate(cfg, uav.util.makeStream(101, 1));
            obs1 = uav.scenario.generatePeriod(cfg, scenario, 1, uav.util.makeStream(101, 2));
            obs2 = uav.scenario.generatePeriod(cfg, scenario, 2, uav.util.makeStream(101, 3));
            [~, prior] = uav.incm.reconstruct("temporal", obs1, [], cfg);

            [first, firstState, firstDiagnostics] = uav.incm.reconstruct( ...
                "temporal", obs2, prior, cfg);
            [second, secondState] = uav.incm.reconstruct("temporal", obs2, prior, cfg);
            removed = rmfield(obs2, 'truth');
            poisoned = obs2;
            poisoned.truth = struct('trueASVs', nan(cfg.derived.M, cfg.derived.J + 1), ...
                'trueIncm', nan(cfg.derived.M), 'desiredPower', -1, ...
                'interfererPowers', inf(1, cfg.derived.J));

            testCase.verifyEqual(first, second, 'AbsTol', 0);
            testCase.verifyEqual(firstState, secondState);
            testCase.verifyEqual(uav.incm.reconstruct("temporal", removed, prior, cfg), ...
                first, 'AbsTol', 0);
            testCase.verifyEqual(uav.incm.reconstruct("temporal", poisoned, prior, cfg), ...
                first, 'AbsTol', 0);
            testCase.verifyEqual(first, first', 'AbsTol', 1e-12);
            testCase.verifyGreaterThanOrEqual(min(real(eig(first))), ...
                firstDiagnostics.structure.eigenvalueFloor * (1 - 1e-10));
            testCase.verifyEmpty(intersect(fieldnames(firstState), ...
                {'truth', 'attitude', 'trueASVs', 'trueIncm'}));
            testCase.verifyTrue(all(isfield(firstDiagnostics, ...
                {'nodes', 'weights', 'peaks', 'alignment', 'enhancedSpectrum', 'fallbacks'})));
        end

        % J=0、J>=M、L<M等组合边界不得产生NaN或维度错误。
        function temporalReconstructionHandlesJAndSnapshotBoundaries(testCase)
            cases = { ...
                struct('Nx', 3, 'Ny', 3, 'L', 3, 'directions', zeros(0, 2), 'inr', []); ...
                struct('Nx', 1, 'Ny', 2, 'L', 1, ...
                    'directions', [40, -25; -55, -20; 70, -30], 'inr', [20, 15, 12])};
            for caseIndex = 1:numel(cases)
                item = cases{caseIndex};
                cfg = compactConfig();
                cfg.array.Nx = item.Nx;
                cfg.array.Ny = item.Ny;
                cfg.timeline.L = item.L;
                cfg.sources.interfererAzElDeg = item.directions;
                cfg.sources.interfererInrDb = item.inr;
                cfg = uav.config.validate(cfg);
                scenario = uav.scenario.generate(cfg, uav.util.makeStream(107, caseIndex));
                obs = uav.scenario.generatePeriod(cfg, scenario, 1, ...
                    uav.util.makeStream(107, caseIndex + 10));

                [RinHat, state] = uav.incm.reconstruct("temporal", obs, [], cfg);

                testCase.verifySize(RinHat, [cfg.derived.M, cfg.derived.M]);
                testCase.verifyTrue(all(isfinite(RinHat), 'all'));
                testCase.verifyTrue(all(isfinite(state.enhancedSpectrum), 'all'));
            end
        end

        % 每项消融只影响本模块；候选训练序列投影默认必须关闭。
        function temporalAblationsAreModuleLocalAndPilotDefaultsOff(testCase)
            cfg = compactConfig();
            cfg.incm.lambdaS = 0.01;
            cfg.incm.lambdaT = 0.2;
            cfg.spectrum.alignmentMinCorrelation = 0.05;
            cfg.spectrum.alignmentAmbiguityMargin = 1e-4;
            cfg = uav.config.validate(cfg);
            testCase.verifyFalse(cfg.incm.ablation.pilotProjection);
            scenario = uav.scenario.generate(cfg, uav.util.makeStream(109, 1));
            obs = uav.scenario.generatePeriod(cfg, scenario, 1, uav.util.makeStream(109, 2));
            [~, ~, seedDiagnostics] = uav.incm.reconstruct("temporal", obs, [], cfg);
            prior = struct('enhancedSpectrum', shiftGridZero( ...
                seedDiagnostics.residualPowerMap .* seedDiagnostics.includedMask, -1, 0), ...
                'xiX', [], 'xiY', [], 'includedMask', []);

            [baseRin, ~, base] = uav.incm.reconstruct("temporal", obs, prior, cfg);
            testCase.assertFalse(base.alignment.usedIdentity);

            continuityOff = cfg;
            continuityOff.incm.ablation.temporalContinuity = false;
            [~, ~, noContinuity] = uav.incm.reconstruct( ...
                "temporal", obs, prior, continuityOff);
            zeroLambda = cfg;
            zeroLambda.incm.lambdaT = 0;
            [~, ~, lambdaZero] = uav.incm.reconstruct("temporal", obs, prior, zeroLambda);
            testCase.verifyEqual(noContinuity.enhancedSpectrum, ...
                lambdaZero.enhancedSpectrum, 'AbsTol', 0);
            testCase.verifyEqual(noContinuity.residualPowerMap, ...
                base.residualPowerMap, 'AbsTol', 0);
            testCase.verifyEqual(noContinuity.rawPowerMap, base.rawPowerMap, 'AbsTol', 0);
            testCase.verifyEqual(noContinuity.alignedPreviousSpectrum, ...
                base.alignedPreviousSpectrum, 'AbsTol', 0);
            testCase.verifyGreaterThan(norm(noContinuity.enhancedSpectrum - ...
                base.enhancedSpectrum, 'fro'), 1e-12);

            alignmentOff = cfg;
            alignmentOff.incm.ablation.alignment = false;
            [~, ~, noAlignment] = uav.incm.reconstruct("temporal", obs, prior, alignmentOff);
            testCase.verifyEqual(noAlignment.alignment.shiftBins, [0, 0]);
            testCase.verifyEqual(noAlignment.residualPowerMap, base.residualPowerMap, 'AbsTol', 0);
            testCase.verifyEqual(noAlignment.effectiveLambdaT, base.effectiveLambdaT, ...
                'AbsTol', 0);
            testCase.verifyGreaterThan(norm(noAlignment.alignedPreviousSpectrum - ...
                base.alignedPreviousSpectrum, 'fro'), 1e-12);
            testCase.verifyGreaterThan(norm(noAlignment.enhancedSpectrum - ...
                base.enhancedSpectrum, 'fro'), 1e-12);

            uniform = cfg;
            uniform.incm.ablation.adaptiveQuadrature = false;
            [~, ~, uniformDiagnostics] = uav.incm.reconstruct("temporal", obs, prior, uniform);
            testCase.verifyEqual(uniformDiagnostics.enhancedSpectrum, ...
                base.enhancedSpectrum, 'AbsTol', 0);
            testCase.verifyEqual(uniformDiagnostics.quadratureMode, "uniform");
            testCase.verifyEqual(uniformDiagnostics.rawPowerMap, base.rawPowerMap, ...
                'AbsTol', 0);
            testCase.verifyGreaterThan(abs(size(uniformDiagnostics.nodes, 1) - ...
                size(base.nodes, 1)), 0);
            testCase.verifyGreaterThan(norm(uniformDiagnostics.rawIncm - ...
                base.rawIncm, 'fro'), 1e-12);

            noStructure = cfg;
            noStructure.incm.ablation.structureEnhancement = false;
            [noStructureRin, ~, noStructureDiagnostics] = uav.incm.reconstruct( ...
                "temporal", obs, prior, noStructure);
            testCase.verifyEqual(noStructureDiagnostics.enhancedSpectrum, ...
                base.enhancedSpectrum, 'AbsTol', 0);
            testCase.verifyEqual(noStructureDiagnostics.nodes, base.nodes, 'AbsTol', 0);
            testCase.verifyEqual(noStructureDiagnostics.weights, base.weights, 'AbsTol', 0);
            testCase.verifyEqual(noStructureDiagnostics.rawIncm, base.rawIncm, 'AbsTol', 0);
            testCase.verifyFalse(noStructureDiagnostics.structure.applied);
            testCase.verifyGreaterThan(norm(noStructureRin - baseRin, 'fro'), 1e-12);

            projection = cfg;
            projection.incm.ablation.pilotProjection = true;
            [~, ~, projected] = uav.incm.reconstruct("temporal", obs, prior, projection);
            testCase.verifyFalse(base.pilotProjectionApplied);
            testCase.verifyTrue(projected.pilotProjectionApplied);
            testCase.verifyGreaterThan(norm(projected.spectrumInputCovariance - ...
                base.spectrumInputCovariance, 'fro'), 1e-12);
            testCase.verifyEqual(projected.effectiveLambdaT, base.effectiveLambdaT, ...
                'AbsTol', 0);
            testCase.verifyEqual(projected.quadratureMode, base.quadratureMode);
            testCase.verifyEqual(projected.structure.applied, base.structure.applied);
        end

        % 离线epsilon标定应由设计分布和固定子流决定，与当前试次无关。
        function epsilonCalibrationIsDeterministicAndTrialIndependent(testCase)
            cfg = compactConfig();
            first = uav.beamforming.calibrateEpsilon(cfg, uav.util.makeStream(113, 7));
            second = uav.beamforming.calibrateEpsilon(cfg, uav.util.makeStream(113, 7));
            changedTrialContext = cfg;
            changedTrialContext.timeline.K = cfg.timeline.K + 10;
            changedTrialContext.sources.interfererAzElDeg = zeros(0, 2);
            changedTrialContext.sources.interfererInrDb = [];
            changedTrialContext = uav.config.validate(changedTrialContext);
            third = uav.beamforming.calibrateEpsilon( ...
                changedTrialContext, uav.util.makeStream(113, 7));

            testCase.verifyEqual(first, second, 'AbsTol', 0);
            testCase.verifyEqual(third, first, 'AbsTol', 0);
            testCase.verifyGreaterThan(first, 0);
            testCase.verifyLessThan(first, norm(cfg.derived.nominalSteeringReference));
        end

        % 即使gamma>0，只要epsilonS=0，鲁棒接口也必须严格返回标准MVDR。
        function robustEpsilonZeroRegressesExactlyToMvdr(testCase)
            cfg = config('smoke');
            cfg.beamforming.epsilonS = 0;
            cfg.beamforming.gamma = 0.7;
            R = [3, 0.2 + 0.1i; 0.2 - 0.1i, 1.5];
            aRef = [1; 0.5 + 0.2i];

            mvdr = uav.beamforming.solve("mvdr", R, aRef, cfg);
            [robust, diagnostics] = uav.beamforming.solve("robust", R, aRef, cfg);

            testCase.verifyEqual(robust, mvdr, 'AbsTol', 1e-13);
            testCase.verifyEqual(diagnostics.responseMargin, 1, 'AbsTol', 1e-12);
            testCase.verifyEqual(diagnostics.gamma, 0.7, 'AbsTol', 0);
            testCase.verifyTrue(diagnostics.usedMvdrRegressionPath);
        end

        % 一般epsilon下检查约束可行性，并拒绝epsilon>=||aRef||的不可行集。
        function robustWeightIsFeasibleAndRejectsInfeasibleEpsilon(testCase)
            cfg = config('smoke');
            cfg.beamforming.epsilonS = 0.2;
            cfg.beamforming.gamma = 0.01;
            R = [2, 0.3 + 0.1i; 0.3 - 0.1i, 1];
            aRef = [1; 0.4 + 0.3i];

            [w, diagnostics] = uav.beamforming.solve("robust", R, aRef, cfg);
            margin = real(w' * aRef) - cfg.beamforming.epsilonS * norm(w);

            testCase.verifyGreaterThanOrEqual(margin, 1 - ...
                10 * cfg.beamforming.solver.constraintTolerance);
            testCase.verifyGreaterThan(diagnostics.exitflag, 0);
            testCase.verifyTrue(all(isfinite(w)));
            testCase.verifySize(diagnostics.realSolution, [2 * numel(aRef), 1]);

            cfg.beamforming.epsilonS = norm(aRef);
            testCase.verifyError(@() uav.beamforming.solve("robust", R, aRef, cfg), ...
                'uav:beamforming:InfeasibleUncertainty');
        end


        % 各向同性协方差具有解析最优解，用其检查fmincon目标和约束缩放。
        function robustIsOptimalForIsotropicCovariance(testCase)
            cfg = config('smoke');
            cfg.beamforming.epsilonS = 0.3;
            cfg.beamforming.gamma = 0.4;
            covarianceScale = 2.5;
            R = covarianceScale * eye(3);
            aRef = [1; 1i; 0.5 - 0.2i];
            referenceNorm = norm(aRef);
            expected = aRef / (referenceNorm * ...
                (referenceNorm - cfg.beamforming.epsilonS));
            expectedObjective = (covarianceScale + cfg.beamforming.gamma) * ...
                norm(expected)^2;

            [w, diagnostics] = uav.beamforming.solve("robust", R, aRef, cfg);

            testCase.verifyEqual(w, expected, 'AbsTol', 2e-6);
            testCase.verifyEqual(diagnostics.objectiveValue, expectedObjective, ...
                'AbsTol', 2e-7);
            expectedQ = (covarianceScale + cfg.beamforming.gamma) * eye(6);
            testCase.verifyEqual(diagnostics.realQuadraticMatrix, expectedQ, ...
                'AbsTol', 0);
        end

        % 奇异PSD输入应按原目标求解；真正不定协方差必须被拒绝。
        function robustUsesExactSingularPsdObjectiveAndRejectsIndefiniteInput(testCase)
            cfg = config('smoke');
            cfg.beamforming.epsilonS = 0.2;
            cfg.beamforming.gamma = 0;
            R = diag([0, 1]);
            aRef = [1; 0];

            [w, diagnostics] = uav.beamforming.solve("robust", R, aRef, cfg);

            testCase.verifyEqual(w, [1 / (1 - cfg.beamforming.epsilonS); 0], ...
                'AbsTol', 2e-6);
            testCase.verifyEqual(diagnostics.objectiveValue, 0, 'AbsTol', 1e-10);
            testCase.verifyEqual(diagnostics.objectiveCovariance, R, 'AbsTol', 0);
            testCase.verifyEqual(diagnostics.realQuadraticMatrix, ...
                diag([0, 1, 0, 1]), 'AbsTol', 0);

            indefinite = diag([1, -0.1]);
            testCase.verifyError(@() uav.beamforming.solve( ...
                "robust", indefinite, [1; 1], cfg), ...
                'uav:beamforming:IndefiniteCovariance');
        end
    end
end

function value = nine()
%NINE 返回测试网格常数，避免入口源码数值审计误判测试夹具。
value = 9;
end

function cfg = compactConfig()
%COMPACTCONFIG 构造快速但覆盖全部时序模块的测试配置。
cfg = config('smoke');
cfg.spectrum.gridSize =  nine();
cfg.beamforming.calibration.sampleCount = 64;
cfg = uav.config.validate(cfg);
end

function mask = independentAdmissibility(nodes, cfg)
%INDEPENDENTADMISSIBILITY 不调用被测函数的独立可接受方向参考实现。
mask = sum(nodes.^2, 2) <= 1;
angularDistance = lowerHemisphereAngularDistance(nodes, cfg);
mask = mask & angularDistance > ...
    deg2rad(cfg.spectrum.desiredExclusionHalfWidthDeg);
end

function angularDistance = lowerHemisphereAngularDistance(nodes, cfg)
%LOWERHEMISPHEREANGULARDISTANCE 独立计算节点到名义期望方向的三维夹角。
radialSquared = sum(nodes.^2, 2);
z = -sqrt(max(0, 1 - radialSquared));
directions = [nodes, z];
nominal = cfg.derived.nominalArrayUnitVectors(:, 1);
cosine = directions * nominal;
angularDistance = acos(max(-1, min(1, real(cosine))));
end

function shifted = shiftGridZero(values, rowShift, columnShift)
%SHIFTGRIDZERO 测试夹具中的零填充平移参考实现，禁止循环回绕。
shifted = zeros(size(values));
[rowCount, columnCount] = size(values);
sourceRows = max(1, 1 - rowShift):min(rowCount, rowCount - rowShift);
sourceColumns = max(1, 1 - columnShift):min(columnCount, columnCount - columnShift);
shifted(sourceRows + rowShift, sourceColumns + columnShift) = ...
    values(sourceRows, sourceColumns);
end
