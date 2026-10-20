classdef TestLinYuePowerModules < matlab.unittest.TestCase
    %TESTLINYUEPOWERMODULES 同候选功率重构、低快拍边界与真值隔离检查。
    methods(TestMethodSetup)
        function paths(testCase)
            root=fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,'src')));
        end
    end
    methods(Test)
        function exactCovarianceRecoversYuePowers(testCase)
            cfg=config('smoke');A=cfg.derived.nominalSteering;p=[10;1000;1000];
            Rx=(A.*p.')*A'+eye(9);
            [R,d]=uav.incm.literaturePointPowers('yue_ls',Rx,A,1,1,0,cfg);
            expected=(A(:,2:end).*p(2:end).')*A(:,2:end)'+eye(9);
            testCase.verifyEqual(d.powers,p,'RelTol',1e-10);
            testCase.verifyEqual(R,expected,'AbsTol',1e-8);
            [N,~]=uav.incm.literaturePointPowers('nnls',Rx,A,1,1,0,cfg);
            testCase.verifyEqual(N,expected,'AbsTol',1e-8);
        end
        function linLiteralDoesNotSilentlyAddNoise(testCase)
            cfg=config('smoke');A=cfg.derived.nominalSteering;
            [R,d]=uav.incm.literaturePointPowers('lin_literal',eye(9),A,1,1,0,cfg);
            testCase.verifyEqual(R,A(:,2:end)*A(:,2:end)'/9,'AbsTol',1e-12);
            testCase.verifyEqual(d.noisePower,0);
            testCase.verifyError(@()uav.incm.literaturePointPowers( ...
                'lin_literal',diag([ones(1,8),0]),A,1,1,0,cfg),'uav:study:NotApplicable');
        end
        function negativeLsPowerIsNotClipped(testCase)
            cfg=config('smoke');A=eye(3);Rx=diag([2,0.5,3]);
            [R,d]=uav.incm.literaturePointPowers('yue_ls',Rx,A,1,1,0,cfg);
            testCase.verifyEqual(d.powers,[1;-0.5;2],'AbsTol',1e-12);
            testCase.verifyEqual(d.negativePowerCount,1);
            testCase.verifyEqual(R,diag([1,0.5,3]),'AbsTol',1e-12);
        end
        function strongerRidgeShrinksFixedDictionaryPowers(testCase)
            cfg=config('smoke');A=eye(3);Rx=diag([11,101,201]);
            [~,a]=uav.incm.literaturePointPowers('nnls',Rx,A,1,1,0,cfg);
            [~,b]=uav.incm.literaturePointPowers('nnls',Rx,A,1,1,10,cfg);
            testCase.verifyEqual(b.powers,a.powers/11,'AbsTol',1e-10);
        end
        function refinedCandidatesRemainTruthFreeAndOutsidePrior(testCase)
            cfg=config('smoke');s=uav.scenario.generate(cfg,uav.util.makeStream(cfg.random.baseSeed,1));
            obs=uav.scenario.generatePeriod(cfg,s,1,uav.util.makeStream(cfg.random.baseSeed,2));
            safe=uav.experiment.algorithmObservation(obs);
            [~,d]=uav.incm.upaCovfit(safe,cfg);
            [A,nodes,r]=uav.spectrum.refineInterfererCandidates(safe,d,cfg);
            testCase.verifyFalse(r.usedTruth);
            testCase.verifyEqual(nodes(1:r.desiredCount,:),d.desiredCandidateNodes);
            testCase.verifySize(A,[9,size(nodes,1)]);
            testCase.verifyTrue(all(uav.spectrum.admissibleDirections(nodes(r.desiredCount+1:end,:),cfg)));
        end
    end
end
