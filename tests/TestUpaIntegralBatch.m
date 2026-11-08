classdef TestUpaIntegralBatch < matlab.unittest.TestCase
    % 验证收敛诊断的代数加速与原始Rec逐点积分严格一致。
    methods(TestMethodSetup)
        function paths(t)
            root=fileparts(fileparts(mfilename('fullpath')));
            t.applyFixture(matlab.unittest.fixtures.PathFixture(root));
            t.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,'src')));
        end
    end
    methods(Test)
        function highPowerSourcesMatchOriginalFormula(t)
            % 高INR窄峰下，检查三角展开没有因消去误差而改变原始积分。
            c=config("smoke");c.array.Nx=3;c.array.Ny=3;c.spectrum.gridSize=101;
            c.signal.snrDb=30;c=uav.config.validate(c);
            A=c.derived.nominalSteering;
            R=(A.*reshape(c.derived.sourcePowers,1,[]))*A'+c.derived.noisePower*eye(c.derived.M);
            settings=struct('nodeChunkSize',509,'matrixBatchSize',2);
            fast=uav.incm.integrateUpaBatch(R,c,101,"endpoint",settings);
            direct=uav.incm.reconstruct("incm_integral_origin",struct('Rx',R),[],c);
            t.verifyLessThan(norm(fast-direct,'fro')/norm(direct,'fro'),1e-10);
        end
        function matchesProductionAndDirectMidpoint(t)
            c=config("smoke");c.array.Nx=3;c.array.Ny=2;c=uav.config.validate(c);
            stream=uav.util.makeStream(c.random.baseSeed,42);
            X=randn(stream,c.derived.M,50)+1i*randn(stream,c.derived.M,50);
            R=X*X'/50;batch=cat(3,R,2*R);
            settings=struct('nodeChunkSize',37,'matrixBatchSize',1);
            for rule=["endpoint","midpoint"]
                G=21;
                if rule=="endpoint",axisValues=linspace(-1,1,G);h=2/(G-1);
                else,h=2/G;axisValues=-1+h/2+(0:G-1)*h;end
                [x,y]=meshgrid(axisValues);nodes=[x(:),y(:)];
                nodes=nodes(uav.spectrum.admissibleDirections(nodes,c),:);
                A=uav.array.steeringUPA(nodes.',3,2,c.array.spacingOverLambda);
                direct=uav.incm.guLeshem2012Rec(R,A,repmat(h^2,1,size(A,2)));
                [fast,d]=uav.incm.integrateUpaBatch(batch,c,G,rule,settings);
                t.verifyLessThan(norm(fast(:,:,1)-direct,'fro')/norm(direct,'fro'),1e-12);
                t.verifyEqual(fast(:,:,2),2*fast(:,:,1),'RelTol',1e-12);
                t.verifyEqual(d.includedNodeCount,size(A,2));
                if rule=="endpoint"
                    origin=uav.incm.reconstruct("incm_integral_origin",struct('Rx',R),[],c);
                    t.verifyEqual(fast(:,:,1),origin,'AbsTol',1e-12);
                end
            end
        end
    end
end
