function [Rh,d]=ulaStatic(obs,thetaReference,cfg,halfWidthDeg)
%ULASTATIC 将旧static的加载、谱底、噪声和结构处理移植到一维角域。
%   不是原式Rec。UPA继续调用既有reconstruct("static")，不改其定义。
g=cfg.experiments.recComprehensive; M=cfg.derived.M;
step=pi/g.ulaGridSize; theta=-pi/2+step/2+step*(0:g.ulaGridSize-1);
A=uav.array.steeringULA(theta,M,cfg.array.spacingOverLambda(2));
scale=max(real(trace(obs.Rx))/M,eps);
loading=cfg.spectrum.diagonalLoading*scale;
loaded=(obs.Rx+obs.Rx')/2+loading*eye(M);
power=1./real(sum(conj(A).*(loaded\A),1));
eigenvalues=sort(real(eig((obs.Rx+obs.Rx')/2)));
zeroTolerance=eps(max(1,max(abs(eigenvalues))))*max(M,size(obs.X,2));
positive=eigenvalues(eigenvalues>zeroTolerance);
floorPower=cfg.incm.noiseFloorRatio*scale;
if size(obs.X,2)<M
    noise=cfg.signal.noisePower;
elseif isempty(positive)
    noise=floorPower;
else
    noise=max(median(positive(1:ceil(cfg.spectrum.noiseEigenvalueFraction*numel(positive)))),floorPower);
end
keep=abs(theta-thetaReference)>deg2rad(halfWidthDeg);
residual=max(power(keep)-(noise+loading)/M,0);
raw=(A(:,keep).*(step*residual))*A(:,keep)'+noise*eye(M);
[Rh,structure]=uav.incm.projectStructure(raw,cfg,true);
d=struct('noisePower',noise,'loading',loading,'structure',structure,'usedTruth',false);
end
