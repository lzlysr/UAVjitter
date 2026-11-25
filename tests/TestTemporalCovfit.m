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
                "incm_covfit_refined_temporal_hybrid","incm_covfit_refined_temporal_capon_base", ...
                "incm_covfit_refined_temporal_pure_subspace"];
            cfg.experiments.classicJitterMethodNames = methods;
            cfg = uav.config.validate(cfg);
            scene = uav.scenario.generate(cfg,uav.util.makeStream(100,1));
            first = uav.scenario.generatePeriod(cfg,scene,1,uav.util.makeStream(100,2));
            second = uav.scenario.generatePeriod(cfg,scene,2,uav.util.makeStream(100,3));

            [firstRecord, state] = uav.experiment.runClassicMethods(first,[],cfg,methods);
            testCase.verifyEqual(firstRecord.weights(:,1:4),repmat(firstRecord.weights(:,1),1,4));
            [secondRecord, next] = uav.experiment.runClassicMethods(second,state,cfg,methods);
            testCase.verifyTrue(secondRecord.diagnostics{2}.reconstruction.temporalApplied);
            reversed = uav.experiment.runClassicMethods(second,state,cfg,fliplr(methods));
            testCase.verifyEqual(secondRecord.weights,fliplr(reversed.weights));
            testCase.verifyNotEmpty(next.temporalCaponCovfit.nodes);
            testCase.verifyNotEmpty(next.temporalHybridCovfit.nodes);
            testCase.verifyNotEmpty(next.temporalSubspaceCaponBase.nodes);
            testCase.verifyNotEmpty(next.temporalPureSubspace.nodes);
            safe = uav.experiment.algorithmObservation(second);
            capon = secondRecord.diagnostics{1}.reconstruction;
            subspace = secondRecord.diagnostics{4}.reconstruction;
            caponNodes = uav.spectrum.fullCaponNodes(safe,capon.noisePower,cfg);
            subspaceNodes = uav.spectrum.shiftInvariantNodes(safe.Rx,subspace.noisePower,cfg);
            testCase.verifyEqual(capon.proposedCandidateNodes,caponNodes);
            testCase.verifyEqual(subspace.proposedCandidateNodes,subspaceNodes);
            testCase.verifyFalse(secondRecord.diagnostics{2}.reconstruction.usedTruth);
        end

        function confirmedTrackCanStayInsideForTwoPeriods(testCase)
            cfg = config("preliminary");
            nominal = cfg.derived.nominalDirectionCosines(:,1).';
            outside = [0.52,0.36];
            inside = [0.64,0.24;0.64,0.20];
            testCase.assertTrue(uav.spectrum.admissibleDirections(outside,cfg));
            testCase.assertFalse(any(uav.spectrum.admissibleDirections(inside,cfg)));
            first = observation(outside);
            [~,state] = uav.incm.temporalCovfit(first,[],cfg);
            testCase.assertTrue(any(state.confirmedOutside));
            second = observation(inside(1,:));
            [~,next,d] = uav.incm.temporalCovfit(second,state,cfg);
            testCase.verifyEqual(d.reclassifiedDesiredPeak,inside(1,:),'AbsTol',1e-12);
            testCase.verifyEqual(d.desiredCandidateNodes,nominal);
            testCase.verifyEqual(next.nodes,inside(1,:),'AbsTol',0.003);
            testCase.verifyTrue(next.confirmedOutside);
            third = observation(inside(2,:));
            [~,last,d] = uav.incm.temporalCovfit(third,next,cfg);
            testCase.verifyEqual(d.reclassifiedDesiredPeak,inside(2,:),'AbsTol',1e-12);
            testCase.verifyEqual(last.nodes,inside(2,:),'AbsTol',0.003);
            testCase.verifyTrue(last.confirmedOutside);

            % 同一冲突缺乏区外确认时，仍保留期望字典的临时峰。
            provisional = state;
            provisional.confirmedOutside(:) = false;
            [~,~,d] = uav.incm.temporalCovfit(second,provisional,cfg);
            testCase.verifyEmpty(d.reclassifiedDesiredPeak);
            testCase.verifyFalse(d.temporalApplied);
            cfg.incm.temporalCovfit.resolveDesiredPeakConflict = false;
            [~,~,legacy] = uav.incm.temporalCovfit(second,state,cfg);
            testCase.verifyEmpty(legacy.reclassifiedDesiredPeak);
            testCase.verifyFalse(legacy.temporalApplied);

            function obs = observation(node)
                % 用理论SCM隔离身份冲突，快拍矩阵仅提供维度。
                A = uav.array.steeringUPA([nominal;node].',cfg.array.Nx,cfg.array.Ny, ...
                    cfg.array.spacingOverLambda);
                M = cfg.derived.M;
                obs = struct('Rx',(A.*[10,400])*A'+eye(M),'X',zeros(M,32));
            end
        end

        function pureSubspaceRunsWithoutGridAndKeepsIdentity(testCase)
            cfg = config("preliminary");
            cfg.derived = rmfield(cfg.derived,'spectrumGrid');
            method = "incm_covfit_refined_temporal_pure_subspace";
            out = [0.52,0.36];
            in = [0.64,0.24];
            nominal = cfg.derived.nominalDirectionCosines(:,1).';
            obs = observation([nominal;out]);
            [baseR,~] = uav.incm.subspaceCovfit(obs,cfg);
            [firstR,state] = uav.incm.reconstruct(method,obs,[],cfg);
            testCase.verifyEqual(firstR,baseR);
            testCase.assertTrue(any(state.confirmedOutside));
            obs = observation([-0.4,0.3;in]);
            [~,next,d] = uav.incm.reconstruct(method,obs,state,cfg);
            testCase.verifyEqual(d.reclassifiedDesiredPeak,in,'AbsTol',1e-10);
            testCase.verifyTrue(any(next.confirmedOutside & vecnorm(next.nodes-in,2,2)<1e-10));

            % 制造无法关联的历史，纯子空间版仍走单周期底座。
            state.nodes = [-0.9,-0.3];
            [baseR,~] = uav.incm.subspaceCovfit(obs,cfg);
            [fallbackR,~,d] = uav.incm.reconstruct(method,obs,state,cfg);
            testCase.verifyEqual(fallbackR,baseR);
            testCase.verifyFalse(d.temporalApplied);

            function obs = observation(nodes)
                A = uav.array.steeringUPA(nodes.',cfg.array.Nx,cfg.array.Ny,cfg.array.spacingOverLambda);
                obs = struct('Rx',(A.*[10,400])*A'+eye(cfg.derived.M),'X',zeros(cfg.derived.M,32));
            end
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
