classdef TestIncmNamesAndSectors < matlab.unittest.TestCase
    %TESTINCMNAMESANDSECTORS 核查方法改名兼容和ULA/UPA独立排除参数。
    methods(TestMethodSetup)
        function addProjectPath(testCase)
            root=fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,'src')));
        end
    end
    methods(Test)
        function historicalNamesNormalizeWithoutChangingOrder(testCase)
            old=["gu_leshem2012_rec","static","CovFit","CovFit_local"];
            names=["incm_integral_origin","incm_integral_processed","incm_covfit","incm_covfit_refined"];
            testCase.verifyEqual(uav.config.canonicalMethodNames(old),names);
            testCase.verifyEqual(uav.config.canonicalMethodNames(old.'),names.');
            cfg=config('smoke');cfg.experiments.classicJitterMethodNames=old;
            cfg.experiments.exp02MethodNames=["oracle","static"];
            cfg=uav.config.validate(cfg);
            testCase.verifyEqual(cfg.experiments.classicJitterMethodNames,names);
            testCase.verifyEqual(cfg.experiments.exp02MethodNames,["oracle","incm_integral_processed"]);
            cfg.experiments.classicSnapshotMethodNames=["CovFit","incm_covfit"];
            testCase.verifyError(@()uav.config.validate(cfg),'uav:config:InvalidExperimentParameter');
        end
        function aliasesProduceSameWeightsInBothArrays(testCase)
            old=["gu_leshem2012_rec","static","CovFit","CovFit_local"];
            names=uav.config.canonicalMethodNames(old);
            for shape={[3,3],[9,1]}
                cfg=config('smoke');cfg.array.Nx=shape{1}(1);cfg.array.Ny=shape{1}(2);
                cfg.timeline.L=64;cfg=uav.config.validate(cfg);
                scene=uav.scenario.generate(cfg,uav.util.makeStream(cfg.random.baseSeed,51));
                obs=uav.scenario.generatePeriod(cfg,scene,1,uav.util.makeStream(cfg.random.baseSeed,52));
                a=uav.experiment.runClassicMethods(obs,[],cfg,old);
                b=uav.experiment.runClassicMethods(obs,[],cfg,names);
                testCase.verifyEqual(a.methodNames,names);
                testCase.verifyEqual(a.weights,b.weights,'AbsTol',1e-12);
                testCase.verifyEqual(a.sinrLinear,b.sinrLinear,'AbsTol',1e-12);
                safe=uav.experiment.algorithmObservation(obs);
                for q=1:4
                    [R,state,d]=uav.incm.reconstruct(old(q),safe,[],cfg);
                    testCase.verifyEqual(state.method,names(q));
                    testCase.verifyEqual(d.method,names(q));
                    testCase.verifyEqual(R,uav.incm.reconstruct(names(q),safe,[],cfg),'AbsTol',1e-12);
                    testCase.verifyEqual(d.desiredExclusionHalfWidthDeg,cfg.derived.desiredExclusionHalfWidthDeg);
                end
            end
        end
        function ulaUsesItsOwnWidthAndReportsGeometry(testCase)
            cfg=config('smoke');cfg.array.Nx=9;cfg.array.Ny=1;cfg=uav.config.validate(cfg);
            testCase.verifyEqual(cfg.spectrum.desiredExclusionHalfWidthDeg,25);
            testCase.verifyEqual(cfg.derived.desiredExclusionHalfWidthDeg,5);
            testCase.verifyEqual(cfg.derived.ulaGeometry.interferersInsideDesiredSector,[false,false]);
            nodes=[cfg.derived.nominalDirectionCosines(1,2:end).',zeros(2,1)];
            testCase.verifyTrue(all(uav.spectrum.admissibleDirections(nodes,cfg)));
            cfg.spectrum.ulaDesiredExclusionHalfWidthDeg=45;cfg=uav.config.validate(cfg);
            testCase.verifyEqual(cfg.derived.ulaGeometry.interferersInsideDesiredSector,[true,false]);
            testCase.verifyEqual(uav.spectrum.admissibleDirections(nodes,cfg),[false;true]);
        end
        function oneDimensionalMaskUsesExactSameHalfWidth(testCase)
            cfg=config('smoke');cfg.array.Nx=9;cfg.array.Ny=1;
            cfg.spectrum.ulaDesiredExclusionHalfWidthDeg=7;cfg=uav.config.validate(cfg);
            reference=asin(cfg.derived.nominalDirectionCosines(1,1));
            theta=reference+deg2rad([-7.1,-6.9,6.9,7.1]);
            testCase.verifyEqual(uav.spectrum.admissibleDirections([sin(theta(:)),zeros(4,1)],cfg), ...
                [true;false;false;true]);
        end
        function changingUlaWidthDoesNotAffectUpaReconstruction(testCase)
            cfg=config('smoke');cfg.array.Nx=3;cfg.array.Ny=3;cfg=uav.config.validate(cfg);
            obs=struct('Rx',diag(1:9),'X',ones(9,32));
            a=uav.incm.reconstruct("incm_integral_processed",obs,[],cfg);
            changed=cfg;changed.spectrum.ulaDesiredExclusionHalfWidthDeg=1;
            changed=uav.config.validate(changed);
            b=uav.incm.reconstruct("incm_integral_processed",obs,[],changed);
            testCase.verifyEqual(a,b,'AbsTol',0);
            testCase.verifyEqual(changed.derived.desiredExclusionHalfWidthDeg,25);
        end
        function historicalConfigurationsKeepTheirOldWidth(testCase)
            cfg=config('smoke');cfg.array.Nx=9;cfg.array.Ny=1;
            % 历史配置也缺少新方向约定标记，因此仍按旧三维投影验证宽度回退。
            cfg.array=rmfield(cfg.array,'ulaDirectionConvention');
            cfg.spectrum=rmfield(cfg.spectrum,'ulaDesiredExclusionHalfWidthDeg');
            cfg=uav.config.validate(cfg);
            testCase.verifyEqual(cfg.derived.desiredExclusionHalfWidthDeg,25);
            testCase.verifyEqual(cfg.derived.ulaGeometry.interferersInsideDesiredSector,[true,true]);
        end
        function invalidUlaWidthsAreRejected(testCase)
            cfg=config('smoke');
            for width=[-1,90,Inf,NaN]
                cfg.spectrum.ulaDesiredExclusionHalfWidthDeg=width;
                testCase.verifyError(@()uav.config.validate(cfg),'uav:config:InvalidAlgorithmParameter');
            end
        end
    end
end
