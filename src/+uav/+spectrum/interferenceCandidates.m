function candidates=interferenceCandidates(obs,cfg)
%INTERFERENCECANDIDATES 两种Capon对照共用观测峰与噪声估计，不拟合源功率。
%   复用现有CovFit的局部极大值/NMS与低特征值噪声估计，以免把前端
%   差异误当成后端收益。只返回干扰候选；期望候选不参与重构或导向校正。
%   settings与其它算法独立；ULA沿用ulaGridSize，UPA沿用gridSize。
cfg=uav.spectrum.useUlaComparisonGrid(cfg);
settings=cfg.incm.caponCandidates;front=cfg;
front.incm.localCovfit.loadingRatio=settings.loadingRatio;
for field={'loadingRatio','peakThresholdRatio','peakNmsRadiusBins','maximumCandidateCount'}
    front.incm.staticCovfit.(field{1})=settings.(field{1});
end
original=uav.incm.covfitFrontEnd(obs,front);
[A,D,nodes,selection]=uav.spectrum.localPeakCandidates(obs.Rx,original,front);
candidates=struct('manifold',A(:,D+1:end),'nodes',nodes(D+1:end,:), ...
    'noisePower',original.noisePower,'noiseSource',original.noiseSource, ...
    'noiseDiagnostics',original.noiseDiagnostics,'selection',selection, ...
    'candidateLoading',original.loading,'candidateGridSize',cfg.spectrum.gridSize, ...
    'usedTruth',false);
end
