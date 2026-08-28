classdef TestTask4TemporalRobust < matlab.unittest.TestCase
    methods (TestMethodSetup)
        function addSourcePath(testCase)
            rootDir = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(rootDir));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(rootDir, 'src')));
        end
    end

    methods (Test)
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
            testCase.verifyGreaterThan(cfg.beamforming.solver.initialInteriorMargin, 0);
            testCase.verifyGreaterThan(cfg.beamforming.solver.feasibilityTolerance, 0);
            testCase.verifyEqual(cfg.beamforming.calibration.quantile, 0.95);
            testCase.verifyGreaterThan(cfg.beamforming.calibration.sampleCount, 0);
        end

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
                'beamforming', 'gamma', -1};
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

        function alignmentUsesIdentityForWeakOrAmbiguousEvidence(testCase)
            cfg = config('smoke');
            weak = ones(7);

            [aligned, diagnostics] = uav.spectrum.align2D(weak, weak, true(7), cfg);

            testCase.verifyEqual(aligned, weak, 'AbsTol', 0);
            testCase.verifyEqual(diagnostics.shiftBins, [0, 0]);
            testCase.verifyTrue(diagnostics.usedIdentity);
            testCase.verifyTrue(contains(diagnostics.fallbackReason, "weak"));
        end

        function temporalProxMatchesHandDerivedBreakpointRegions(testCase)
            p = [0.1; 1; 2.1; 4];
            g = [2; 2; 2; 2];
            lambdaS = [1; 0.4; 0.4; 0.4];
            lambdaT = [0.2; 0.2; 0.8; 0.8];

            c = uav.incm.temporalProx(p, g, lambdaS, lambdaT);

            testCase.verifyEqual(c, [0; 0.9; 2; 3.4], 'AbsTol', 1e-12);
        end

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

        function adaptiveQuadratureHasNoOverlapAndConservesGridArea(testCase)
            cfg = config('smoke');
            cfg.spectrum.gridSize = 11;
            cfg = uav.config.validate(cfg);
            [~, spectrumDiagnostics] = uav.spectrum.loadedCapon(eye(cfg.derived.M), cfg);
            includedMask = spectrumDiagnostics.visibleMask;
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
            testCase.verifyTrue(all(sum(quadrature.nodes.^2, 2) <= 1.25));
        end

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

        function temporalAblationsAreModuleLocalAndPilotDefaultsOff(testCase)
            cfg = compactConfig();
            testCase.verifyFalse(cfg.incm.ablation.pilotProjection);
            scenario = uav.scenario.generate(cfg, uav.util.makeStream(109, 1));
            obs = uav.scenario.generatePeriod(cfg, scenario, 1, uav.util.makeStream(109, 2));
            prior = struct('enhancedSpectrum', zeros(cfg.spectrum.gridSize), ...
                'xiX', [], 'xiY', [], 'includedMask', []);
            prior.enhancedSpectrum(4, 4) = 50;

            [~, ~, base] = uav.incm.reconstruct("temporal", obs, prior, cfg);

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

            alignmentOff = cfg;
            alignmentOff.incm.ablation.alignment = false;
            [~, ~, noAlignment] = uav.incm.reconstruct("temporal", obs, prior, alignmentOff);
            testCase.verifyEqual(noAlignment.alignment.shiftBins, [0, 0]);
            testCase.verifyEqual(noAlignment.residualPowerMap, base.residualPowerMap, 'AbsTol', 0);

            uniform = cfg;
            uniform.incm.ablation.adaptiveQuadrature = false;
            [~, ~, uniformDiagnostics] = uav.incm.reconstruct("temporal", obs, prior, uniform);
            testCase.verifyEqual(uniformDiagnostics.enhancedSpectrum, ...
                base.enhancedSpectrum, 'AbsTol', 0);
            testCase.verifyEqual(uniformDiagnostics.quadratureMode, "uniform");

            noStructure = cfg;
            noStructure.incm.ablation.structureEnhancement = false;
            [~, ~, noStructureDiagnostics] = uav.incm.reconstruct( ...
                "temporal", obs, prior, noStructure);
            testCase.verifyEqual(noStructureDiagnostics.enhancedSpectrum, ...
                base.enhancedSpectrum, 'AbsTol', 0);
            testCase.verifyFalse(noStructureDiagnostics.structure.applied);

            projection = cfg;
            projection.incm.ablation.pilotProjection = true;
            [~, ~, projected] = uav.incm.reconstruct("temporal", obs, prior, projection);
            testCase.verifyFalse(base.pilotProjectionApplied);
            testCase.verifyTrue(projected.pilotProjectionApplied);
            testCase.verifyGreaterThan(norm(projected.spectrumInputCovariance - ...
                base.spectrumInputCovariance, 'fro'), 1e-12);
        end

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

        function robustEpsilonZeroRegressesExactlyToMvdr(testCase)
            cfg = config('smoke');
            cfg.beamforming.epsilonS = 0;
            cfg.beamforming.gamma = 0;
            R = [3, 0.2 + 0.1i; 0.2 - 0.1i, 1.5];
            aRef = [1; 0.5 + 0.2i];

            mvdr = uav.beamforming.solve("mvdr", R, aRef, cfg);
            [robust, diagnostics] = uav.beamforming.solve("robust", R, aRef, cfg);

            testCase.verifyEqual(robust, mvdr, 'AbsTol', 1e-13);
            testCase.verifyEqual(diagnostics.responseMargin, 1, 'AbsTol', 1e-12);
            testCase.verifyTrue(diagnostics.usedMvdrRegressionPath);
        end

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
    end
end

function value = nine()
value = 9;
end

function cfg = compactConfig()
cfg = config('smoke');
cfg.spectrum.gridSize =  nine();
cfg.beamforming.calibration.sampleCount = 64;
cfg = uav.config.validate(cfg);
end
