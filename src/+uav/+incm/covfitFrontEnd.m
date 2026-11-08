function original = covfitFrontEnd(obs,cfg)
%COVFITFRONTEND 两种CovFit共用的期望候选、加载谱与噪声估计。
%   保持原细化版的期望字典和噪声算法，不读取真实方向或真实干扰数。
%   干扰局部峰由后续共用localPeakCandidates提取；本函数不提前拟合功率。
if ~isstruct(obs) || isfield(obs,'truth') || ~isfield(obs,'X') || ...
        ~isfield(obs,'Rx') || ~ismatrix(obs.X) || isempty(obs.X) || ...
        any(~isfinite(obs.X),'all')
    error('uav:incm:InvalidObservation','CovFit必须使用剥离真值后的有限快拍和协方差。');
end
M=cfg.derived.M;
if size(obs.X,1)~=M || ~isequal(size(obs.Rx),[M,M]) || any(~isfinite(obs.Rx),'all')
    error('uav:incm:InvalidObservation','CovFit快拍和协方差维数必须与阵元数一致。');
end
local=cfg;
local.spectrum.diagonalLoading=cfg.incm.localCovfit.loadingRatio;
[powerMap,spectrum]=uav.spectrum.loadedCapon(obs.Rx,local);
allNodes=[spectrum.xiX(:),spectrum.xiY(:)];
outside=uav.spectrum.admissibleDirections(allNodes,cfg);
desiredMask=spectrum.visibleMask & ~reshape(outside,size(powerMap));
nominal=cfg.derived.nominalDirectionCosines(:,1).';
if cfg.derived.arrayType=="ULA",nominal(2)=0;end
desired=nominal;
indices=find(desiredMask & isfinite(powerMap));
if isempty(indices)
    desiredFallback="nominal-only";
else
    [~,order]=max(real(powerMap(indices)));index=indices(order);
    node=[spectrum.xiX(index),spectrum.xiY(index)];
    if all(max(abs(desired-node),[],2)>cfg.spectrum.quadratureDedupTolerance)
        desired(end+1,:)=node;
        desiredFallback="strongest-finite";
    else
        desiredFallback="deduplicated-with-nominal";
    end
end
snapshotCount=size(obs.X,2);
if snapshotCount<M
    noisePower=cfg.signal.noisePower;
    noiseSource="configured-underdetermined";
    noiseDiagnostics=struct('eigenvalues',zeros(0,1),'lowerEigenvalueCount',0);
else
    % 与原细化版一致：较小特征值的中位数估噪，不改成真值噪声。
    eigenvalues=sort(real(eig((obs.Rx+obs.Rx')/2)),'ascend');
    lowerCount=max(1,ceil(M*cfg.spectrum.noiseEigenvalueFraction));
    rawEstimate=median(eigenvalues(1:lowerCount));
    scale=max(real(trace(obs.Rx))/M,eps);
    noisePower=max(rawEstimate,eps(scale));noiseSource="lower-eigenvalues";
    noiseDiagnostics=struct('eigenvalues',eigenvalues,'lowerEigenvalueCount',lowerCount);
end
original=struct('desiredCandidateNodes',desired,'noisePower',noisePower, ...
    'noiseSource',noiseSource,'noiseDiagnostics',noiseDiagnostics, ...
    'effectiveSnapshotCount',snapshotCount,'loading',spectrum.loading, ...
    'desiredCandidateFallback',desiredFallback,'usedTruth',false);
end
