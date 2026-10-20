classdef TestRecComprehensive < matlab.unittest.TestCase
    %TESTRECCOMPREHENSIVE 检验独立几何、配对复现和不适用边界。
    methods(Test)
        function planIsDeterministicAndRich(testCase)
            cfg=config('preliminary'); [a,c]=uav.experiment.recComprehensivePlan(cfg);
            [b,d]=uav.experiment.recComprehensivePlan(cfg);
            testCase.verifyEqual(a,b);testCase.verifyEqual(c,d);
            testCase.verifyEqual(sum([a.randomIndex]>0),64);
            testCase.verifyEqual(sum(string({a.family})=="ula_cartesian"),200);
            testCase.verifyEqual(sum(string({a.family})=="upa_cartesian"),90);
            testCase.verifyEqual(numel(cfg.experiments.recComprehensive.snapshotCounts),11);
            testCase.verifyEqual(sum(string({c.block})=="snapshots"),128);
        end
        function zeroJitterAndSnapshotStatuses(testCase)
            cfg=config('smoke');[a,c]=uav.experiment.recComprehensivePlan(cfg);
            for array=["ULA","UPA"]
                id=find(string({c.block})=="snapshots" & string({c.arrayType})==array,1);
                condition=c(id);condition.trials=2;condition.snapshotCounts=[8,32];
                r=uav.experiment.runRecComprehensiveCase(cfg,a(condition.geometryIndex),condition);
                testCase.verifyEqual(r.status(:,1,1),zeros(2,1,'int8'));
                testCase.verifyEqual(r.status(:,2,1),ones(2,1,'int8'));
                testCase.verifyEqual(r.status(:,:,2:end),ones(2,2,3,'int8'));
                testCase.verifyGreaterThanOrEqual(min(r.oracleLossDb,[],'all','omitnan'),-1e-6);
                testCase.verifyEqual(r.desiredGain(r.status==1),ones(sum(r.status(:)==1),1),'AbsTol',1e-7);
                repeated=uav.experiment.runRecComprehensiveCase(cfg,a(condition.geometryIndex),condition);
                testCase.verifyEqual(r.sinrLinear,repeated.sinrLinear);
                testCase.verifyEqual(r.snapshotSubstreams,repeated.snapshotSubstreams);
            end
        end
        function sourceCountIncludesEmptyAndOverloaded(testCase)
            cfg=config('smoke');[a,c]=uav.experiment.recComprehensivePlan(cfg);
            for J=[0,9]
                id=find(string({c.block})=="source_count" & string({c.arrayType})=="UPA" & ...
                    arrayfun(@(q)numel(q.inrDb)==J,c),1);
                r=uav.experiment.runRecComprehensiveCase(cfg,a(c(id).geometryIndex),c(id));
                testCase.verifyEqual(r.status,ones(1,1,4,'int8'));
                testCase.verifyTrue(all(isfinite(r.sinrLinear),'all'));
            end
        end
        function medianIntervalCoverageContract(testCase)
            testCase.verifyEqual(uav.experiment.recMedianInterval(1:50,.95),[18,33]);
            testCase.verifyTrue(all(isnan(uav.experiment.recMedianInterval(1:3,.95))));
            testCase.verifyEqual(uav.experiment.recMedianInterval(ones(1,50),.95),[1,1]);
        end
        function oneDimensionalStaticMatchesDirectIntegral(testCase)
            cfg=config('smoke');cfg.array.Nx=1;cfg.array.Ny=9;cfg=uav.config.validate(cfg);
            obs=struct('X',zeros(9,8),'Rx',eye(9)); reference=0;half=25;
            count=cfg.experiments.recComprehensive.ulaGridSize;step=pi/count;
            theta=-pi/2+step/2+step*(0:count-1);keep=abs(theta)>deg2rad(half);
            B=uav.array.steeringULA(theta(keep),9,cfg.array.spacingOverLambda(2));
            % 白噪声谱恰等于所减底座，积分余量为舍入量级，最后应返回单位阵。
            [r,d]=uav.incm.ulaStatic(obs,reference,cfg,half);
            testCase.verifyEqual(r,eye(9),'AbsTol',1e-12);
            testCase.verifyFalse(d.usedTruth);testCase.verifySize(B,[9,sum(keep)]);
        end
    end
end
