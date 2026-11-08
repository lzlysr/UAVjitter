classdef TestUpaJacobianIntegral < matlab.unittest.TestCase
    % 在研究分支验证测度换元、同节点权重和直接积分，不涉及方向估计。
    methods(TestMethodSetup)
        function paths(t)
            root=fileparts(fileparts(mfilename('fullpath')));
            t.applyFixture(matlab.unittest.fixtures.PathFixture(root));
            t.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,'src')));
        end
    end
    methods(Test)
        function azimuthElevationUsesConsistentDensity(t)
            % 等方位/高度角与等方位/z坐标必须给出相同立体角积分。
            c=config("smoke");c.array.Nx=3;c.array.Ny=2;
            c.spectrum.desiredExclusionHalfWidthDeg=0;c=uav.config.validate(c);
            s=c.experiments.upaGridConvergence;G=401;M=c.derived.M;
            [R,d]=uav.incm.integrateUpaBatch(eye(M),c,G,"azel_midpoint",s, ...
                ["dxi_x_dxi_y","solid_angle","angle_product"]);
            reference=uav.incm.integrateUpaBatch(eye(M),c,G,"surface_midpoint",s, ...
                ["dxi_x_dxi_y","solid_angle"]);
            t.verifyEqual(d(1).includedArea,pi,'AbsTol',1e-4);
            t.verifyEqual(d(2).includedArea,2*pi,'AbsTol',1e-4);
            t.verifyEqual(d(3).includedArea,pi^2,'AbsTol',1e-10);
            for m=1:2
                t.verifyLessThan(norm(R(:,:,1,m)-reference(:,:,1,m),'fro')/norm(reference(:,:,1,m),'fro'),3e-5);
            end
            % 保持同一个Capon数值却换成角度乘积测度，会改变估计器而非仅换坐标。
            t.verifyGreaterThan(norm(R(:,:,1,3)/trace(R(:,:,1,3))-R(:,:,1,2)/trace(R(:,:,1,2)),'fro'),.05);
        end
        function fullHemisphereAnalyticIntegral(t)
            c=config("smoke");c.array.Nx=3;c.array.Ny=2;
            c.spectrum.desiredExclusionHalfWidthDeg=0;c=uav.config.validate(c);
            s=c.experiments.upaGridConvergence;G=401;M=c.derived.M;
            [R,d]=uav.incm.integrateUpaBatch(eye(M),c,G,"surface_midpoint",s, ...
                ["dxi_x_dxi_y","solid_angle"]);
            % 零排除半宽下，圆盘面积pi；下半球立体角2pi。
            t.verifyEqual(d(1).includedArea,pi,'AbsTol',1e-10);
            t.verifyEqual(d(2).includedArea,2*pi,'AbsTol',1e-10);
            x=kron((0:c.array.Nx-1)',ones(c.array.Ny,1));
            y=repmat((0:c.array.Ny-1)',c.array.Nx,1);
            distance=sqrt(((x-x.')*c.array.spacingOverLambda(1)).^2+ ...
                ((y-y.')*c.array.spacingOverLambda(2)).^2);
            % 对单位白噪声，Capon谱是1/M；半球平面阵响应积分为2pi*sinc(2r)。
            expected=2*pi/M*sinc(2*distance);
            t.verifyLessThan(norm(R(:,:,1,2)-expected,'fro')/norm(expected,'fro'),2e-5);
        end
        function sameNodesMatchDirectOriginalKernel(t)
            c=config("smoke");c.array.Nx=3;c.array.Ny=3;c=uav.config.validate(c);
            s=c.experiments.upaGridConvergence;G=51;M=c.derived.M;
            A=c.derived.nominalSteering;R=A*diag(c.derived.sourcePowers)*A'+eye(M);
            for rule=["surface_midpoint","surface_shifted"]
                [fast,d]=uav.incm.integrateUpaBatch(R,c,G,rule,s,["dxi_x_dxi_y","solid_angle"]);
                phi=-pi+pi/G+(0:G-1)*2*pi/G;
                if rule=="surface_shifted",phi=phi-pi/G;end
                z=-1+1/(2*G)+(0:G-1)/G;[az,z]=meshgrid(phi,z);
                radius=sqrt(1-z(:).^2);nodes=[radius.*cos(az(:)),radius.*sin(az(:))];
                keep=uav.spectrum.admissibleDirections(nodes,c);nodes=nodes(keep,:);z=z(keep);
                D=uav.array.steeringUPA(nodes.',3,3,c.array.spacingOverLambda);
                weights=[abs(z),ones(size(z))]*2*pi/G^2;
                for m=1:2
                    direct=uav.incm.guLeshem2012Rec(R,D,weights(:,m));
                    t.verifyLessThan(norm(fast(:,:,1,m)-direct,'fro')/norm(direct,'fro'),1e-10);
                    t.verifyEqual(d(m).includedNodeCount,size(D,2));
                end
            end
        end
    end
end
