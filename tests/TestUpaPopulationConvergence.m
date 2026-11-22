classdef TestUpaPopulationConvergence < matlab.unittest.TestCase
    % 原式一致性、半球Jacobian解析控制与独立产物管理。
    methods (TestMethodSetup)
        function paths(t)
            root=fileparts(fileparts(mfilename('fullpath')));
            for p={root,fullfile(root,'src'),fullfile(root,'experiments')}
                t.applyFixture(matlab.unittest.fixtures.PathFixture(p{1}));
            end
        end
    end
    methods (Test)
        function optionalNativeMatchesMatlab(t)
            % 无编译器的机器可正常使用纯MATLAB路径；启用MEX时逐元素核对。
            t.assumeEqual(exist('upaPopulationQuadratureNative','file'),3);
            c=config("smoke");c.array.Nx=9;c.array.Ny=9;c.signal.snrDb=30;
            c=uav.config.validate(c);A=c.derived.nominalSteering;
            R=A*diag(c.derived.sourcePowers)*A'+eye(c.derived.M);
            inputs=cat(3,R,eye(c.derived.M));s=c.experiments.upaPopulationConvergence;
            for rule=["endpoint","surface_midpoint","surface_shifted","azel_midpoint"]
                if rule=="endpoint",measures="dxi_x_dxi_y";
                else,measures=["dxi_x_dxi_y","solid_angle"];end
                for G=[51,401]
                    s.nativeAcceleration=false;
                    [a,da]=uav.incm.integrateUpaBatch(inputs,c,G,rule,s,measures);
                    s.nativeAcceleration=true;
                    [b,db]=uav.incm.integrateUpaBatch(inputs,c,G,rule,s,measures);
                    t.verifyLessThan(norm(a(:)-b(:))/norm(a(:)),s.parityTolerance);
                    t.verifyEqual([da.includedNodeCount],[db.includedNodeCount]);
                    t.verifyEqual([da.includedArea],[db.includedArea],'AbsTol',1e-9);
                end
            end
        end
        function nineByNineMatchesProductionFormula(t)
            c=config("smoke");c.array.Nx=9;c.array.Ny=9;c.spectrum.gridSize=31;
            c.signal.snrDb=30;c=uav.config.validate(c);A=c.derived.nominalSteering;
            R=A*diag(c.derived.sourcePowers)*A'+c.signal.noisePower*eye(c.derived.M);
            s=c.experiments.upaPopulationConvergence;s.nodeChunkSize=113;
            fast=uav.incm.integrateUpaBatch(R,c,31,"endpoint",s);
            direct=uav.incm.reconstruct("incm_integral_origin",struct('Rx',R),[],c);
            t.verifyLessThan(norm(fast-direct,'fro')/norm(direct,'fro'),s.parityTolerance);
        end
        function sphericalAndAngularNodesMatchDirectSum(t)
            c=config("smoke");c.array.Nx=3;c.array.Ny=2;c=uav.config.validate(c);
            G=21;s=c.experiments.upaPopulationConvergence;s.nodeChunkSize=37;
            A=c.derived.nominalSteering;R=A*diag(c.derived.sourcePowers)*A'+eye(c.derived.M);
            for rule=["surface_midpoint","surface_shifted","azel_midpoint"]
                phi=-pi+pi/G+(0:G-1)*2*pi/G;
                if rule=="surface_shifted",phi=phi-pi/G;end
                if rule=="azel_midpoint"
                    alpha=-pi/2+pi/(4*G)+(0:G-1)*pi/(2*G);
                    [az,el]=meshgrid(phi,alpha);rr=cos(el(:));z=sin(el(:));
                    weights=[abs(z).*rr,rr]*pi^2/G^2;
                else
                    [az,z]=meshgrid(phi,-1+1/(2*G)+(0:G-1)/G);
                    z=z(:);rr=sqrt(1-z.^2);weights=[abs(z),ones(G^2,1)]*2*pi/G^2;
                end
                nodes=[rr.*cos(az(:)),rr.*sin(az(:))];
                keep=uav.spectrum.admissibleDirections(nodes,c);
                nodes=nodes(keep,:);weights=weights(keep,:);
                D=uav.array.steeringUPA(nodes.',3,2,c.array.spacingOverLambda);
                fast=uav.incm.integrateUpaBatch(R,c,G,rule,s,["dxi_x_dxi_y","solid_angle"]);
                for m=1:2
                    direct=uav.incm.guLeshem2012Rec(R,D,weights(:,m));
                    t.verifyLessThan(norm(fast(:,:,1,m)-direct,'fro')/norm(direct,'fro'),1e-9);
                end
            end
        end
        function whiteHemisphereHasAnalyticMassAndCovariance(t)
            c=config("smoke");c.array.Nx=3;c.array.Ny=2;
            c.spectrum.desiredExclusionHalfWidthDeg=0;c=uav.config.validate(c);
            s=c.experiments.upaPopulationConvergence;M=c.derived.M;
            [R,d]=uav.incm.integrateUpaBatch(eye(M),c,401,"surface_midpoint",s, ...
                ["dxi_x_dxi_y","solid_angle"]);
            t.verifyEqual(d(1).includedArea,pi,'AbsTol',1e-10);
            t.verifyEqual(d(2).includedArea,2*pi,'AbsTol',1e-10);
            x=kron((0:2)',ones(2,1));y=repmat((0:1)',3,1);
            r=sqrt(((x-x')*c.array.spacingOverLambda(1)).^2+ ...
                ((y-y')*c.array.spacingOverLambda(2)).^2);
            expected=2*pi/M*sinc(2*r);
            t.verifyLessThan(norm(R(:,:,1,2)-expected,'fro')/norm(expected,'fro'),3e-5);
        end
        function standaloneRunnerDoesNotOverwriteEarlierRun(t)
            folder=t.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            c=config("smoke");c.output.rootDir=folder.Folder;
            s=c.experiments.upaPopulationConvergence;s.arrayShape=[3,2];s.snrDb=[-10,10];
            s.gridSizes=[11,21,41];s.extensionGridSizes=81;s.parityGridSize=11;
            c.experiments.upaPopulationConvergence=s;
            a=diagnose_upa_population_convergence(c);
            % 续算也必须另建目录；复用已完成节点，不重新生成原结果。
            c.experiments.upaPopulationConvergence.resumeCheckpoint= ...
                fullfile(a.outputDirectory,'checkpoint.mat');
            b=diagnose_upa_population_convergence(c);
            t.verifyNotEqual(a.outputDirectory,b.outputDirectory);
            t.verifyTrue(isfile(fullfile(a.outputDirectory,'study.mat')));
            t.verifyTrue(isfile(fullfile(b.outputDirectory,'study.mat')));
            t.verifyEqual(a.cfg.attitude.stdDeg,[0,0,0]);
            t.verifyEqual(a.metadata.actualCounts.monteCarloTrials,0);
            t.verifyTrue(all(a.parity.relativeMatrixError<s.parityTolerance));
            t.verifyEqual(a.records(1).R,b.records(1).R,'AbsTol',1e-12);
            t.verifyEqual(numel(a.records),numel(b.records));
        end
    end
end
