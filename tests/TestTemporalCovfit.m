classdef TestTemporalCovfit < matlab.unittest.TestCase
    % 验证时序CovFit的初始化、观测候选和配对统计。

    methods (TestMethodSetup)
        function addPaths(testCase)
            root = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,'src')));
        end
    end

    methods (Test)
        function firstPeriodMatchesBaselineAndStateResets(testCase)
            cfg = config("smoke");
            scenario = uav.scenario.generate(cfg,uav.util.makeStream(100,1));
            first = uav.scenario.generatePeriod(cfg,scenario,1,uav.util.makeStream(100,2));
            second = uav.scenario.generatePeriod(cfg,scenario,2,uav.util.makeStream(100,3));
            methods = ["incm_covfit_refined","incm_covfit_refined_temporal"];
            [firstRecord, state] = uav.experiment.runClassicMethods(first,[],cfg,methods);
            testCase.verifyEqual(firstRecord.weights(:,1),firstRecord.weights(:,2));
            testCase.verifyFalse(firstRecord.diagnostics{2}.reconstruction.temporalApplied);
            testCase.verifyNotEmpty(state.temporalCovfit.nodes);
            [nextRecord, next] = uav.experiment.runClassicMethods(second,state,cfg,methods);
            reset = uav.experiment.runClassicMethods(second,[],cfg,methods);
            testCase.verifyEqual(reset.weights(:,1),reset.weights(:,2));
            diagnostics = nextRecord.diagnostics{2}.reconstruction;
            testCase.verifyTrue(diagnostics.temporalApplied);

            % 直接移除导频字段，时序重构及其递推状态应逐元素相同。
            safe = uav.experiment.algorithmObservation(second);
            [expectedR, expectedState] = uav.incm.temporalCovfit(safe,state.temporalCovfit,cfg);
            [actualR, actualState] = uav.incm.temporalCovfit( ...
                rmfield(safe,'pilot'),state.temporalCovfit,cfg);
            testCase.verifyEqual(actualR,expectedR);
            testCase.verifyEqual(actualState,expectedState);
            testCase.verifyFalse(nextRecord.diagnostics{2}.reconstruction.usedTruth);
            testCase.verifyGreaterThan(min(nextRecord.incmMinimumEigenvalue),0);
            testCase.verifyNotEmpty(next.temporalCovfit.powers);
        end

        function fullCaponRetainsPeakInsideDesiredExclusion(testCase)
            cfg = config("smoke");
            insideNode = cfg.derived.nominalDirectionCosines(:,1).';
            nodes = [insideNode;-0.4,0.3];
            A = uav.array.steeringUPA(nodes.',cfg.array.Nx,cfg.array.Ny, ...
                cfg.array.spacingOverLambda);
            Rx = (A.*[1000,600])*A' + eye(cfg.derived.M);
            obs = struct('Rx',Rx);
            [full,powers] = uav.spectrum.fullCaponNodes(obs,1,cfg);
            mask = uav.spectrum.admissibleDirections(full,cfg);
            testCase.verifyTrue(any(~mask & vecnorm(full-insideNode,2,2)<0.002));
            testCase.verifyTrue(all(isfinite(powers) & powers>=0));
        end

        function hybridCaponUsesTheDesignatedFrontendInEachZone(testCase)
            cfg = config("smoke");
            inside = cfg.derived.nominalDirectionCosines(:,1).';
            outside = [-0.4,0.3];
            A = uav.array.steeringUPA([inside;outside].',cfg.array.Nx,cfg.array.Ny, ...
                cfg.array.spacingOverLambda);
            Rx = (A.*[1000,600])*A' + eye(cfg.derived.M);
            obs = struct('Rx',Rx);
            subspace = uav.spectrum.shiftInvariantNodes(Rx,1,cfg);
            capon = uav.spectrum.fullCaponNodes(obs,1,cfg);
            expectedInside = subspace(~uav.spectrum.admissibleDirections(subspace,cfg),:);
            expectedOutside = capon(uav.spectrum.admissibleDirections(capon,cfg),:);
            testCase.verifyNotEmpty(expectedInside);
            testCase.verifyNotEmpty(expectedOutside);
            [hybrid,powers] = uav.spectrum.hybridCaponNodes(obs,1,cfg);
            testCase.verifyEqual(hybrid,[expectedInside;expectedOutside]);
            testCase.verifyTrue(all(isfinite(powers) & powers>0));
        end

        function frontendsKeepIndependentHistory(testCase)
            cfg = config("smoke");
            methods = ["incm_covfit_refined_temporal","incm_covfit_refined_temporal_capon", ...
                "incm_covfit_refined_temporal_hybrid"];
            cfg.experiments.classicJitterMethodNames = methods;
            cfg = uav.config.validate(cfg);
            scene = uav.scenario.generate(cfg,uav.util.makeStream(100,1));
            first = uav.scenario.generatePeriod(cfg,scene,1,uav.util.makeStream(100,2));
            second = uav.scenario.generatePeriod(cfg,scene,2,uav.util.makeStream(100,3));

            [firstRecord, state] = uav.experiment.runClassicMethods(first,[],cfg,methods);
            testCase.verifyEqual(firstRecord.weights,repmat(firstRecord.weights(:,1),1,3));
            [secondRecord, next] = uav.experiment.runClassicMethods(second,state,cfg,methods);
            testCase.verifyTrue(secondRecord.diagnostics{2}.reconstruction.temporalApplied);
            reversed = uav.experiment.runClassicMethods(second,state,cfg,fliplr(methods));
            testCase.verifyEqual(secondRecord.weights,fliplr(reversed.weights));
            testCase.verifyNotEmpty(next.temporalCaponCovfit.nodes);
            testCase.verifyNotEmpty(next.temporalHybridCovfit.nodes);
            testCase.verifyFalse(secondRecord.diagnostics{2}.reconstruction.usedTruth);
        end

        function planarSubspaceRecoversKnownModes(testCase)
            cfg = config("smoke");
            nodes = [0.4,0.3;0.2,-0.5;-0.3,0.3];
            powers = [10;200;500];
            noisePower = 1;
            A = uav.array.steeringUPA(nodes.',cfg.array.Nx,cfg.array.Ny, ...
                cfg.array.spacingOverLambda);
            Rx = (A.*powers.')*A' + noisePower*eye(cfg.derived.M);
            [estimated, estimatedPowers] = uav.spectrum.shiftInvariantNodes(Rx,noisePower,cfg);
            [estimatedPowers, order] = sort(estimatedPowers);
            testCase.verifyEqual(estimated(order,:),nodes,'AbsTol',1e-10);
            testCase.verifyEqual(estimatedPowers,powers,'AbsTol',1e-8);
        end

        function temporalPowerPenaltyMovesTowardPreviousPower(testCase)
            cfg = config("smoke");
            A = cfg.derived.nominalSteering;
            a = A(:,2);
            noise = 1;
            Rx = 5*(a*a') + noise*eye(cfg.derived.M);
            penalty = cfg.derived.M^2;
            q = uav.incm.fitNonnegativePowers(Rx,a,noise,0,penalty,3);
            testCase.verifyEqual(q,4,'AbsTol',1e-12);
            kkt = uav.incm.nonnegativeFitKkt(Rx,a,noise,q,0,penalty,3);
            testCase.verifyLessThan(kkt,1e-12);
        end

        function pairedIntervalsPreserveTrialPairing(testCase)
            cfg = config("smoke");
            samples = reshape(1:18,3,3,2);
            samples(:,:,2) = 2*samples(:,:,1);
            result = uav.experiment.summarizeTrials(samples,cfg,1,[1,2]);
            testCase.verifyEqual(result.pairedDifferenceDb,10*log10(2),'AbsTol',1e-12);
            testCase.verifyEqual(result.pairedDifferenceCi95Db,repmat(10*log10(2),2,1), ...
                'AbsTol',1e-12);
        end
    end
end
