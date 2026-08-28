classdef TestTask1Contracts < matlab.unittest.TestCase
    methods (TestMethodSetup)
        function addSourcePath(testCase)
            rootDir = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(rootDir));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(...
                fullfile(rootDir, 'src')));
        end
    end

    methods (Test)
        function configSuppliesApprovedDefaultsAndDerivedQuantities(testCase)
            cfg = config('smoke');

            testCase.verifyEqual(cfg.array.Nx, 3);
            testCase.verifyEqual(cfg.array.Ny, 3);
            testCase.verifyEqual(cfg.array.spacingOverLambda, [0.5, 0.5]);
            testCase.verifyEqual(cfg.sources.desiredAzElDeg, [0, -35]);
            testCase.verifyEqual(cfg.sources.interfererAzElDeg, [40, -25; -55, -20]);
            testCase.verifyEqual(cfg.sources.interfererInrDb, [20, 15]);
            testCase.verifyEqual(cfg.signal.snrDb, 10);
            testCase.verifyEqual(cfg.timeline.dtSec, 1e-3);
            testCase.verifyEqual(cfg.timeline.L, 8);
            testCase.verifyEqual(cfg.attitude.rho, [0.95, 0.95, 0.98]);
            testCase.verifyEqual(cfg.random.baseSeed, 20260827);
            testCase.verifyEqual(cfg.experiments.jitterSweepDeg, [0, 0.5, 1, 2, 4, 8]);
            testCase.verifyEqual(cfg.experiments.snapshotSweep, [2, 4, 8, 16, 32, 64]);
            testCase.verifyEqual(cfg.derived.M, 9);
            testCase.verifyEqual(cfg.derived.J, 2);
            testCase.verifyEqual(cfg.derived.noisePower, 1);
            testCase.verifyEqual(cfg.derived.desiredPower, 10, 'AbsTol', 1e-12);
            testCase.verifyEqual(cfg.derived.interfererPowers, [100, 10^(1.5)], ...
                'AbsTol', 1e-12);
            testCase.verifyEqual(cfg.derived.worldUnitVectors(:, 1), ...
                [0.819152044288992; 0; -0.573576436351046], 'AbsTol', 1e-12);
            testCase.verifySize(cfg.derived.nominalSteeringReference, [9, 1]);
            testCase.verifyEqual(norm(cfg.derived.nominalSteeringReference), 3, ...
                'AbsTol', 1e-12);
            testCase.verifyTrue(isfolder(cfg.output.cacheDir));
            testCase.verifyTrue(isfolder(cfg.output.rawDir));
        end

        function profilesOnlyChangeWorkloadAndGridResolution(testCase)
            smoke = config('smoke');
            preliminary = config('preliminary');
            paper = config('paper');

            testCase.verifyEqual([smoke.monteCarlo.nTrials, preliminary.monteCarlo.nTrials, ...
                paper.monteCarlo.nTrials], [3, 100, 500]);
            testCase.verifyLessThan(smoke.timeline.K, preliminary.timeline.K);
            testCase.verifyLessThan(preliminary.timeline.K, paper.timeline.K);
            testCase.verifyLessThan(smoke.spectrum.gridSize, preliminary.spectrum.gridSize);
            testCase.verifyLessThan(preliminary.spectrum.gridSize, paper.spectrum.gridSize);

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

        function configRejectsUnknownProfile(testCase)
            testCase.verifyError(@() config('draft'), 'uav:config:UnknownProfile');
        end

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

        function runAllSupportsAQuietOrchestrationModeForTests(testCase)
            result = run_all('smoke', 'RunTests', false, 'RunExperiments', false);

            testCase.verifyEqual(result.profile, "smoke");
            testCase.verifyEmpty(result.tests);
            testCase.verifyEmpty(result.experiments);
        end
    end
end
