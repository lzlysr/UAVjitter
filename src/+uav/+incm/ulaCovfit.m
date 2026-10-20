function [Rh,d] = ulaCovfit(obs,thetaReference,cfg,halfWidthDeg,refine)
%ULACOVFIT 将Capon候选加非负点源拟合思路移植到一维ULA角域。
%   本方法不是原论文Rec，也不是直接调用UPA实现；无真实方向或干扰数输入。
%   thetaReference单位弧度，只用于分区与期望干扰项的联合拟合。
g = cfg.experiments.recUlaGeometry;
M = g.sensorCount;
step = pi/g.candidateGridSize;
theta = -pi/2+step/2+step*(0:g.candidateGridSize-1);
A = uav.array.steeringULA(theta,M,g.spacingOverLambda);
loaded = (obs.Rx+obs.Rx')/2 + g.loadingRatio*real(trace(obs.Rx))/M*eye(M);
power = 1./real(sum(conj(A).*(loaded\A),1));
inside = abs(theta-thetaReference)<=deg2rad(halfWidthDeg);
[~,relative] = max(power(inside)); indices = find(inside);
desired = [thetaReference,theta(indices(relative))];
% 复用原有“阈值＋网格NMS”语义，候选个数不使用真实干扰数。
peakCfg = cfg;
peakCfg.derived.M = M;
peakCfg.spectrum.peakThresholdRatio = g.peakThresholdRatio;
peakCfg.spectrum.peakNmsRadiusBins = g.peakNmsRadiusBins;
maps = struct('xiX',zeros(size(theta)),'xiY',sin(theta),'visibleMask',true(size(theta)));
[~,peaks] = uav.spectrum.selectPeaks(power,maps,~inside,peakCfg);
interferers = theta(peaks.linearIndices(1:min(numel(peaks.linearIndices),g.maximumCandidateCount)));
if refine
    options = optimoptions('fmincon','Algorithm','sqp','Display','off', ...
        'SpecifyObjectiveGradient',true, ...
        'MaxIterations',cfg.experiments.recRootCause.candidateRefinement.maxIterations, ...
        'OptimalityTolerance',cfg.experiments.recRootCause.candidateRefinement.optimalityTolerance, ...
        'StepTolerance',cfg.experiments.recRootCause.candidateRefinement.stepTolerance);
    inverse = loaded\eye(M);
    desired(2) = refineOne(desired(2),true);
    for k = 1:numel(interferers)
        interferers(k) = refineOne(interferers(k),false);
    end
end
desired = dedup(desired,[],g.dedupToleranceRad);
interferers = dedup(interferers,desired,g.dedupToleranceRad);
nodes = [desired,interferers];
B = uav.array.steeringULA(nodes,M,g.spacingOverLambda);
if size(obs.X,2)<M
    noisePower = cfg.signal.noisePower; % 欠采样时沿用允许的配置噪声先验。
else
    eigenvalues = sort(real(eig((obs.Rx+obs.Rx')/2)));
    noisePower = max(median(eigenvalues(1:ceil(M*cfg.spectrum.noiseEigenvalueFraction))),eps);
end
[powers,fit] = uav.incm.fitNonnegativePowers(obs.Rx,B,noisePower,g.fitRegularization);
Bi = B(:,numel(desired)+1:end); interferencePowers = powers(numel(desired)+1:end);
Rh = (Bi.*interferencePowers.')*Bi' + noisePower*eye(M);
Rh = (Rh+Rh')/2;
d = struct('nodesRad',nodes,'desiredCount',numel(desired),'powers',powers, ...
    'noisePower',noisePower,'fit',fit,'usedTruth',false);
    function result = refineOne(seed,isDesired)
        half = deg2rad(halfWidthDeg);
        if isDesired
            lower = max(-pi/2,thetaReference-half); upper = min(pi/2,thetaReference+half);
        elseif seed<thetaReference
            lower = -pi/2; upper = thetaReference-half;
        else
            lower = thetaReference+half; upper = pi/2;
        end
        [candidate,value,flag] = fmincon(@objective,seed,[],[],[],[],lower,upper,[],options);
        result = seed;
        if flag>0 && isfinite(value) && value<=objective(seed)
            result = candidate;
        end
    end
    function [value,gradient] = objective(angle)
        a = uav.array.steeringULA(angle,M,g.spacingOverLambda);
        derivative = 1i*2*pi*g.spacingOverLambda*(0:M-1)'*cos(angle).*a;
        value = real(a'*inverse*a);
        gradient = 2*real(derivative'*inverse*a);
    end
end

function output = dedup(nodes,other,tolerance)
%DEDUP 保留先出现的候选，禁止重复指纹造成不必要的拟合病态。
output = zeros(1,0);
for node = nodes(:).'
    if all(abs(node-[other,output])>tolerance), output(end+1) = node; end %#ok<AGROW>
end
end
