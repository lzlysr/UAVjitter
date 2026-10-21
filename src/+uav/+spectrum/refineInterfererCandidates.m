function [A,nodes,d] = refineInterfererCandidates(obs, original, cfg)
%REFINEINTERFERERCANDIDATES 观测驱动局部Capon峰的亚网格细化，不校正期望参考。
%   低加载与峰细化分别作为显式消融，不冒充Lin的BCS或Taylor原算法。
%   仅优化期望区外的干扰节点；候选数量不读取真实J，期望字典保持不变。
[~,D,nodes] = uav.spectrum.localPeakCandidates(obs.Rx,original,cfg);
M=cfg.derived.M;
Rc=(obs.Rx+obs.Rx')/2+cfg.incm.staticCovfit.loadingRatio*real(trace(obs.Rx))/M*eye(M);
inverse=Rc\eye(M);
sensorX=kron((0:cfg.array.Nx-1)',ones(cfg.array.Ny,1));
sensorY=repmat((0:cfg.array.Ny-1)',cfg.array.Nx,1);
phaseDerivative=2*pi*[cfg.array.spacingOverLambda(1)*sensorX,cfg.array.spacingOverLambda(2)*sensorY];
radius=cfg.incm.localCovfit.radiusBins*2/(cfg.spectrum.gridSize-1);
reference=cfg.derived.nominalArrayUnitVectors(:,1);
cosineLimit=cosd(cfg.spectrum.desiredExclusionHalfWidthDeg);
% 主线参数集中在config.m，不依赖研究分支的实验专用字段。
settings=cfg.incm.localCovfit.solver;
options=optimoptions('fmincon','Algorithm','sqp','Display','off', ...
    'SpecifyObjectiveGradient',true,'MaxIterations',settings.maxIterations, ...
    'OptimalityTolerance',settings.optimalityTolerance, ...
    'StepTolerance',settings.stepTolerance, ...
    'ConstraintTolerance',settings.constraintTolerance);
accepted=0;
exitflags=zeros(size(nodes,1)-D,1);acceptedMask=false(size(exitflags));
for q=D+1:size(nodes,1)
    seed=nodes(q,:).';
    [value,score,flag]=fmincon(@objective,seed,[],[],[],[], ...
        max(seed-radius,-1),min(seed+radius,1),@constraints,options);
    c=constraints(value);
    exitflags(q-D)=flag;
    if flag>0 && all(c<=settings.feasibilityTolerance) && score<=objective(seed)
        nodes(q,:)=value.';accepted=accepted+1;
        acceptedMask(q-D)=true;
    end
end
A=uav.array.steeringUPA(nodes.',cfg.array.Nx,cfg.array.Ny,cfg.array.spacingOverLambda);
d=struct('desiredCount',D,'acceptedCount',accepted,'usedTruth',false, ...
    'exitflags',exitflags,'accepted',acceptedMask,'solver',settings);
    function [value,gradient]=objective(xi)
        a=exp(1i*phaseDerivative*xi);
        z=inverse*a;
        value=real(a'*z);
        da=1i*phaseDerivative.*a;
        gradient=2*real(da'*z);
    end
    function [c,ceq]=constraints(xi)
        lowerZ=-sqrt(max(0,1-sum(xi.^2)));
        c=[sum(xi.^2)-1;reference'*[xi;lowerZ]-cosineLimit];ceq=[];
    end
end
