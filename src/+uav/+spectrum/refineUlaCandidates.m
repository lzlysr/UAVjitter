function [A,nodes,diagnostics] = refineUlaCandidates(obs,nodes,desiredCount,cfg)
%REFINEULACANDIDATES 仅在一维等效角内细化干扰峰，不校正期望参考或字典。
%   借鉴研究分支ulaCovfit的连续角搜索；搜索半径和求解器读main配置，
%   不复制研究分支的全分区搜索或期望候选细化。无真实方向/干扰数输入。
M = cfg.derived.M;
Rc = (obs.Rx+obs.Rx')/2 + cfg.incm.staticCovfit.loadingRatio*real(trace(obs.Rx))/M*eye(M);
inverse = Rc\eye(M);
reference = asin(max(-1,min(1,cfg.derived.nominalDirectionCosines(1,1))));
halfWidth = deg2rad(uav.spectrum.desiredExclusionHalfWidth(cfg));
radius = cfg.incm.localCovfit.radiusBins*cfg.derived.spectrumGrid.step;
% radiusBins不自动改参；四方法统一高密度网格后，细化物理半径随步长缩小。
settings = cfg.incm.localCovfit.solver;
options = optimoptions('fmincon','Algorithm','sqp','Display','off', ...
    'SpecifyObjectiveGradient',true,'MaxIterations',settings.maxIterations, ...
    'OptimalityTolerance',settings.optimalityTolerance, ...
    'StepTolerance',settings.stepTolerance,'ConstraintTolerance',settings.constraintTolerance);
exitflags = zeros(size(nodes,1)-desiredCount,1);
accepted = false(size(exitflags));
for q = desiredCount+1:size(nodes,1)
    seed = asin(max(-1,min(1,nodes(q,1))));
    lower = max(-pi/2,seed-radius); upper = min(pi/2,seed+radius);
    if seed < reference
        upper = min(upper,reference-halfWidth);
    else
        lower = max(lower,reference+halfWidth);
    end
    [angle,score,flag] = fmincon(@objective,seed,[],[],[],[],lower,upper,[],options);
    exitflags(q-desiredCount) = flag;
    % 严格保持期望排除规则；边界解不偷偷移动进期望组。
    if flag>0 && isfinite(score) && score<=objective(seed) && ...
            angle>=lower-settings.feasibilityTolerance && angle<=upper+settings.feasibilityTolerance && ...
            uav.spectrum.admissibleDirections([sin(angle),0],cfg)
        nodes(q,:) = [sin(angle),0]; accepted(q-desiredCount) = true;
    end
end
A = uav.array.steeringUPA(nodes.',cfg.array.Nx,cfg.array.Ny,cfg.array.spacingOverLambda);
diagnostics = struct('desiredCount',desiredCount,'acceptedCount',nnz(accepted), ...
    'usedTruth',false,'exitflags',exitflags,'accepted',accepted,'solver',settings, ...
    'coordinate',"theta=asin(u_x)",'radiusRad',radius);
    function [value,gradient] = objective(angle)
        phase = 2*pi*cfg.array.spacingOverLambda(1)*(0:M-1)';
        a = exp(1i*phase*sin(angle)); z = inverse*a;
        derivative = 1i*phase*cos(angle).*a;
        value = real(a'*z); gradient = 2*real(derivative'*z);
    end
end
