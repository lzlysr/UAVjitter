classdef TestTask3BaselinesIncm < matlab.unittest.TestCase
    methods (TestMethodSetup)
        function addSourcePath(testCase)
            rootDir = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(rootDir));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(rootDir, 'src')));
        end
    end

    methods (Test)
        function configSuppliesStaticIncmNumericalParameters(testCase)
            cfg = config('smoke');

            testCase.verifyEqual(cfg.spectrum.noiseEigenvalueFraction, 0.5);
            testCase.verifyEqual(cfg.incm.noiseFloorRatio, 1e-6);
        end

        function validationRejectsInvalidStaticIncmNumericalParameters(testCase)
            cfg = config('smoke');
            cfg.spectrum.noiseEigenvalueFraction = 0;
            testCase.verifyError(@() uav.config.validate(cfg), ...
                'uav:config:InvalidAlgorithmParameter');

            cfg = config('smoke');
            cfg.incm.noiseFloorRatio = -1;
            testCase.verifyError(@() uav.config.validate(cfg), ...
                'uav:config:InvalidAlgorithmParameter');
        end

        function mvdrHasClosedFormDistortionlessResponse(testCase)
            cfg = config('smoke');
            R = diag([2, 3]);
            aRef = [1; 1];

            [w, diagnostics] = uav.beamforming.solve("mvdr", R, aRef, cfg);

            testCase.verifyEqual(w, [3; 2] / 5, 'AbsTol', 1e-12);
            testCase.verifyEqual(w' * aRef, 1, 'AbsTol', 1e-12);
            testCase.verifyEqual(diagnostics.loading, 0, 'AbsTol', 0);
        end

        function dlmvdrProducesFiniteDistortionlessWeightsForSingularScm(testCase)
            cfg = config('smoke');
            R = [1, 1; 1, 1];
            aRef = [1; -1];

            [w, diagnostics] = uav.beamforming.solve("dlmvdr", R, aRef, cfg);

            testCase.verifyTrue(all(isfinite(w)));
            testCase.verifyEqual(w' * aRef, 1, 'AbsTol', 1e-12);
            testCase.verifyGreaterThan(diagnostics.loading, 0);
        end

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
