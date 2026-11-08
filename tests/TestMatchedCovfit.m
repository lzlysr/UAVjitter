classdef TestMatchedCovfit < matlab.unittest.TestCase
    %TESTMATCHEDCOVFIT 验证两种CovFit只有连续方向细化开关不同。
    methods(TestMethodSetup)
        function addPaths(testCase)
            root=fileparts(fileparts(mfilename('fullpath')));
            for path={root,fullfile(root,'src')}
                testCase.applyFixture(matlab.unittest.fixtures.PathFixture(path{1}));
            end
        end
    end
    methods(Test)
        function commonParametersAreBasedOnExistingRefinedVersion(testCase)
            c=config("smoke");
            testCase.verifyEqual(c.incm.localCovfit.loadingRatio,1e-3);
            testCase.verifyEqual(c.incm.localCovfit.fitRegularization,0);
            fields=fieldnames(c.incm.staticCovfit);
            for i=1:numel(fields)
                testCase.verifyEqual(c.incm.staticCovfit.(fields{i}),c.incm.localCovfit.(fields{i}));
            end
            % 旧字段只保留兼容镜像，不能让两方法分别使用不一致的参数。
            c.incm.staticCovfit.loadingRatio=0.12;c.incm.staticCovfit.fitRegularization=5000;
            c.incm.localCovfit.loadingRatio=0.002;c=uav.config.validate(c);
            testCase.verifyEqual(c.incm.staticCovfit.loadingRatio,0.002);
            testCase.verifyEqual(c.incm.staticCovfit.fitRegularization,0);
        end
        function bothArraysShareCandidatesNoiseAndFit(testCase)
            for shape={[3,3],[9,1]}
                for L=[8,32]
                    c=config("smoke");c.array.Nx=shape{1}(1);c.array.Ny=shape{1}(2);
                    c.timeline.L=L;c=uav.config.validate(c);
                    scene=uav.scenario.generate(c,uav.util.makeStream(c.random.baseSeed,211));
                    obs=uav.scenario.generatePeriod(c,scene,1,uav.util.makeStream(c.random.baseSeed,212));
                    safe=uav.experiment.algorithmObservation(obs);
                    [plain,a]=uav.incm.upaCovfit(safe,c);
                    [refined,b]=uav.incm.localCovfit(safe,c);
                    testCase.verifyEqual(a.commonParameters,b.commonParameters);
                    testCase.verifyEqual(a.initialCandidateNodes,b.initialCandidateNodes);
                    testCase.verifyEqual(a.initialCandidateManifold,b.initialCandidateManifold);
                    testCase.verifyEqual(a.noisePower,b.noisePower);
                    testCase.verifyEqual(a.noiseDiagnostics,b.noiseDiagnostics);
                    testCase.verifyEqual(a.loading,b.loading);
                    testCase.verifyEqual(a.peakSelection,b.peakSelection);
                    testCase.verifyEqual(a.desiredCandidateNodes,b.desiredCandidateNodes);
                    testCase.verifyFalse(a.refinementApplied);testCase.verifyTrue(b.refinementApplied);
                    testCase.verifyFalse(a.refinement.enabled);testCase.verifyTrue(b.refinement.enabled);
                    testCase.verifyEqual(a.candidateManifold,a.initialCandidateManifold,'AbsTol',0);
                    testCase.verifyEqual(plain,uav.incm.localCovfit(safe,c,false),'AbsTol',0);
                    testCase.verifyLessThan(a.relativeKktResidual,c.incm.localCovfit.kktTolerance);
                    testCase.verifyLessThan(b.relativeKktResidual,c.incm.localCovfit.kktTolerance);
                    testCase.verifyTrue(all(isfinite(plain),'all') && all(isfinite(refined),'all'));
                    % 独立重做关闭细化分支的NNLS及矩阵合成，排除隐藏的额外处理。
                    p=uav.incm.fitNonnegativePowers(obs.Rx,a.initialCandidateManifold,a.noisePower,c.incm.localCovfit.fitRegularization);
                    D=size(a.desiredCandidateNodes,1);A=a.initialCandidateManifold(:,D+1:end);
                    expected=(A.*p(D+1:end).')*A'+a.noisePower*eye(c.derived.M);
                    testCase.verifyEqual(plain,(expected+expected')/2,'AbsTol',1e-12);
                    record=uav.experiment.runClassicMethods(obs,[],c,["incm_covfit","incm_covfit_refined"]);
                    testCase.verifyEqual(record.snapshotSignatures(1,:),record.snapshotSignatures(2,:));
                    for i=1:2
                        testCase.verifyEqual(record.diagnostics{i}.beamforming.method,"mvdr");
                        testCase.verifyLessThan(abs(record.weights(:,i)'*c.derived.nominalSteeringReference-1),1e-10);
                    end
                end
            end
        end
    end
end
