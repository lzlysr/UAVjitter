classdef TestTask2GeometryScenario < matlab.unittest.TestCase
    methods (TestMethodSetup)
        function addSourcePath(testCase)
            rootDir = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(rootDir));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(rootDir, 'src')));
        end
    end

    methods (Test)
        function rotationZYXIsProperOrthogonal(testCase)
            R = uav.geometry.rotationZYX([pi / 6; -pi / 5; pi / 4]);

            testCase.verifyEqual(R' * R, eye(3), 'AbsTol', 1e-12);
            testCase.verifyEqual(det(R), 1, 'AbsTol', 1e-12);
        end

        function zeroAttitudeKeepsWorldDirectionsAndCoordinates(testCase)
            uWorld = [1, 0; 0, 1; 0, 0];
            uArray = uav.geometry.worldToArray(eye(3), uWorld);
            [xi, azElRad] = uav.geometry.directionCoordinates(uArray);

            testCase.verifyEqual(uArray, uWorld, 'AbsTol', 1e-12);
            testCase.verifyEqual(xi, [1, 0; 0, 1], 'AbsTol', 1e-12);
            testCase.verifyEqual(azElRad, [0, pi / 2; 0, 0], 'AbsTol', 1e-12);
        end

        function nonzeroRotationAndCoordinateConventionMatchesLiteralZYXCase(testCase)
            etaRad = [pi / 6; -pi / 4; pi / 3];
            R = uav.geometry.rotationZYX(etaRad);
            uWorld = [0; 0; -1];
            uArray = uav.geometry.worldToArray(R, uWorld);
            [xi, azElRad] = uav.geometry.directionCoordinates(uArray);
            expectedR = [ ...
                0.353553390593274, -0.926776695296637, 0.126826484044322; ...
                0.612372435695795, 0.126826484044322, -0.780330085889911; ...
                0.707106781186548, 0.353553390593274, 0.612372435695795];

            testCase.verifyEqual(R, expectedR, 'AbsTol', 1e-12);
            testCase.verifyEqual(uArray, [-0.707106781186548; ...
                -0.353553390593274; -0.612372435695795], 'AbsTol', 1e-12);
            testCase.verifyEqual(xi, [-0.707106781186548; -0.353553390593274], ...
                'AbsTol', 1e-12);
            testCase.verifyEqual(azElRad, [-2.677945044588987; -0.659058035826409], ...
                'AbsTol', 1e-12);
        end

        function geometryAndSteeringRejectInvalidPhysicalInputs(testCase)
            testCase.verifyError(@() uav.geometry.worldToArray(diag([1, 1, 2]), [0; 0; -1]), ...
                'uav:geometry:InvalidRotation');
            testCase.verifyError(@() uav.geometry.worldToArray(eye(3), [0; 0; -2]), ...
                'uav:geometry:InvalidDirectionNorm');
            testCase.verifyError(@() uav.geometry.directionCoordinates([1.01; 0; 0]), ...
                'uav:geometry:InvalidDirectionNorm');
            testCase.verifyError(@() uav.array.steeringUPA([1.01; 0], 2, 2, [0.5, 0.5]), ...
                'uav:array:OutsideVisibleDisk');

            [xi, azElRad] = uav.geometry.directionCoordinates([0; 0; 1]);
            a = uav.array.steeringUPA(xi, 2, 2, [0.5, 0.5]);
            testCase.verifyEqual(azElRad, [0; pi / 2], 'AbsTol', 1e-12);
            testCase.verifyEqual(a, ones(4, 1), 'AbsTol', 1e-12);
        end

        function rectangularSteeringUsesVxKronVyOrderAndHasExpectedNorm(testCase)
            a = uav.array.steeringUPA([0.25; -0.5], 2, 3, [0.5, 0.5]);
            expected = [1; -1i; -1; ...
                (1 + 1i) / sqrt(2); (1 - 1i) / sqrt(2); (-1 - 1i) / sqrt(2)];

            testCase.verifyEqual(a, expected, 'AbsTol', 1e-12);
            testCase.verifyEqual(abs(a), ones(6, 1), 'AbsTol', 1e-12);
            testCase.verifyEqual(norm(a), sqrt(6), 'AbsTol', 1e-12);
        end

        function scenarioHasStationaryAttitudeCovarianceAndZeroJitterBehavior(testCase)
            cfg = config('smoke');
            cfg.array.Nx = 1;
            cfg.array.Ny = 1;
            cfg.sources.interfererAzElDeg = zeros(0, 2);
            cfg.sources.interfererInrDb = [];
            cfg.timeline.K = 20000;
            cfg = uav.config.validate(cfg);

            scenario = uav.scenario.generate(cfg, uav.util.makeStream(13, 1));
            empiricalVariance = var(scenario.attitudeOffsetRad, 0, 2);
            targetVariance = diag(scenario.Qeta);

            testCase.verifyEqual(scenario.Qnu, scenario.Qeta - scenario.F * ...
                scenario.Qeta * scenario.F', 'AbsTol', 1e-14);
            testCase.verifyEqual(empiricalVariance, targetVariance, 'RelTol', 0.35);

            cfg.attitude.stdDeg = [0, 0, 0];
            cfg = uav.config.validate(cfg);
            zeroScenario = uav.scenario.generate(cfg, uav.util.makeStream(13, 1));
            testCase.verifyEqual(zeroScenario.attitudeOffsetRad, zeros(3, cfg.timeline.K), ...
                'AbsTol', 0);
            testCase.verifyEqual(zeroScenario.rotations(:, :, 1), cfg.derived.nominalRotation, ...
                'AbsTol', 1e-12);
            testCase.verifyEqual(zeroScenario.steeringTruth(:, :, 1), cfg.derived.nominalSteering, ...
                'AbsTol', 1e-12);
        end

        function periodIsDeterministicAndSeparatesTruthFromAlgorithmInputs(testCase)
            cfg = config('smoke');
            scenarioA = uav.scenario.generate(cfg, uav.util.makeStream(22, 1));
            scenarioB = uav.scenario.generate(cfg, uav.util.makeStream(22, 1));
            obsA = uav.scenario.generatePeriod(cfg, scenarioA, 3, uav.util.makeStream(22, 2));
            obsB = uav.scenario.generatePeriod(cfg, scenarioB, 3, uav.util.makeStream(22, 2));

            testCase.verifyEqual(obsA.X, obsB.X, 'AbsTol', 0);
            testCase.verifyEqual(obsA.Rx, obsB.Rx, 'AbsTol', 0);
            testCase.verifyEqual(obsA.pilot, obsB.pilot, 'AbsTol', 0);
            testCase.verifySize(obsA.X, [cfg.derived.M, cfg.timeline.L]);
            testCase.verifySize(obsA.Rx, [cfg.derived.M, cfg.derived.M]);
            testCase.verifyEqual(obsA.Rx, obsA.X * obsA.X' / cfg.timeline.L, 'AbsTol', 1e-12);
            testCase.verifyFalse(isfield(obsA, 'trueASVs'));
            testCase.verifyFalse(isfield(obsA, 'trueIncm'));
            testCase.verifyTrue(isfield(obsA, 'truth'));
            testCase.verifySize(obsA.truth.trueASVs, [cfg.derived.M, cfg.derived.J + 1]);
            testCase.verifySize(obsA.truth.trueIncm, [cfg.derived.M, cfg.derived.M]);
            testCase.verifyTrue(all(isfinite(obsA.X), 'all'));
            testCase.verifyTrue(all(isfinite(obsA.Rx), 'all'));
        end

        function periodSupportsFewerSnapshotsThanSensorsAndZeroInterferers(testCase)
            cfg = config('smoke');
            cfg.sources.interfererAzElDeg = zeros(0, 2);
            cfg.sources.interfererInrDb = [];
            cfg.timeline.L = 3;
            cfg = uav.config.validate(cfg);
            scenario = uav.scenario.generate(cfg, uav.util.makeStream(29, 1));
            obs = uav.scenario.generatePeriod(cfg, scenario, 1, uav.util.makeStream(29, 2));

            testCase.verifyLessThan(cfg.timeline.L, cfg.derived.M);
            testCase.verifySize(obs.X, [cfg.derived.M, 3]);
            testCase.verifySize(obs.truth.trueASVs, [cfg.derived.M, 1]);
            testCase.verifyEqual(obs.truth.trueIncm, cfg.derived.noisePower * eye(cfg.derived.M), ...
                'AbsTol', 1e-12);
        end

        function periodSupportsMoreInterferersThanSensors(testCase)
            cfg = config('smoke');
            cfg.array.Nx = 1;
            cfg.array.Ny = 2;
            cfg.sources.interfererAzElDeg = [40, -25; -55, -20; 70, -30];
            cfg.sources.interfererInrDb = [20, 15, 12];
            cfg = uav.config.validate(cfg);
            scenario = uav.scenario.generate(cfg, uav.util.makeStream(31, 1));
            obs = uav.scenario.generatePeriod(cfg, scenario, 1, uav.util.makeStream(31, 2));

            testCase.verifyGreaterThanOrEqual(cfg.derived.J, cfg.derived.M);
            testCase.verifySize(obs.truth.trueASVs, [2, 4]);
            testCase.verifyTrue(all(isfinite(obs.truth.trueIncm), 'all'));
        end

        function largeSnapshotScmConvergesToTheoreticalCovariance(testCase)
            cfg = config('smoke');
            cfg.timeline.L = 50000;
            cfg = uav.config.validate(cfg);
            scenario = uav.scenario.generate(cfg, uav.util.makeStream(37, 1));
            obs = uav.scenario.generatePeriod(cfg, scenario, 1, uav.util.makeStream(37, 2));
            A = obs.truth.trueASVs;
            expectedRx = A * diag(obs.truth.sourcePowers) * A' + ...
                cfg.derived.noisePower * eye(cfg.derived.M);

            testCase.verifyEqual(obs.Rx, expectedRx, 'RelTol', 0.08);
        end

        function qpskSignalsAndDerivedNoiseHaveExpectedStatistics(testCase)
            cfg = config('smoke');
            cfg.timeline.L = 50000;
            cfg = uav.config.validate(cfg);
            scenario = uav.scenario.generate(cfg, uav.util.makeStream(43, 1));
            obs = uav.scenario.generatePeriod(cfg, scenario, 1, uav.util.makeStream(43, 2));
            A = obs.truth.trueASVs;
            desired = sqrt(obs.truth.desiredPower) * A(:, 1) * obs.pilot;
            interference = A(:, 2:end) * (sqrt(obs.truth.interfererPowers(:)) .* ...
                obs.truth.interfererSymbols);
            noise = obs.X - desired - interference;
            pilotInterfererCorrelation = mean(obs.pilot .* conj( ...
                obs.truth.interfererSymbols(1, :)));

            testCase.verifyEqual(abs(obs.pilot), ones(1, cfg.timeline.L), 'AbsTol', 1e-12);
            testCase.verifyEqual(abs(obs.truth.interfererSymbols), ...
                ones(cfg.derived.J, cfg.timeline.L), 'AbsTol', 1e-12);
            testCase.verifyLessThan(abs(pilotInterfererCorrelation), 0.02);
            testCase.verifyEqual(mean(abs(noise).^2, 'all'), cfg.derived.noisePower, ...
                'RelTol', 0.03);
            testCase.verifyLessThan(abs(mean(noise.^2, 'all')), 0.01);
        end
    end
end
