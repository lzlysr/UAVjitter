function [RinHat,diagnostics] = localCovfit(obs,cfg,applyRefinement)
%LOCALCOVFIT Capon局部干扰峰细化与非负协方差拟合，固定最终期望参考。
%   移植自研究分支1477a46的refined_nnls0，不是Le全分区版本。
%   流程：低加载谱→局部峰→有界连续干扰方向细化→联合NNLS→剔除期望组。
%   两种CovFit共用本实现和incm.localCovfit参数；仅applyRefinement不同。
%   默认true保持原细化版；false保留相同局部峰及NNLS，仅跳过连续角优化。
%   只细化干扰节点，期望拟合字典不变；不读取任何真实方向或真实INCM。
if nargin<3,applyRefinement=true;end
% 在整个候选/细化链路中统一网格，不能只加密前端却仍用旧步长细化。
cfg=uav.spectrum.useUlaComparisonGrid(cfg);
settings=cfg.incm.localCovfit;
front=cfg;
% 共用细化版参数作为唯一来源；旧staticCovfit字段只供底层接口兼容。
front.incm.staticCovfit.loadingRatio=settings.loadingRatio;
front.incm.staticCovfit.peakThresholdRatio=settings.peakThresholdRatio;
front.incm.staticCovfit.peakNmsRadiusBins=settings.peakNmsRadiusBins;
front.incm.staticCovfit.maximumCandidateCount=settings.maximumCandidateCount;
front.incm.staticCovfit.fitRegularization=settings.fitRegularization;
original=uav.incm.covfitFrontEnd(obs,front);
[A,D,nodes,candidateSelection]=uav.spectrum.localPeakCandidates(obs.Rx,original,front);
initialNodes=nodes;initialManifold=A;
if applyRefinement
    [A,nodes,refinement]=uav.spectrum.refineInterfererCandidates(obs,nodes,D,front);
else
    % 未细化分支不调用连续优化器；后面的拟合、剔除与合成完全相同。
    refinement=struct('desiredCount',D,'acceptedCount',0,'usedTruth',false);
end
refinement.enabled=logical(applyRefinement);
[powers,fit]=uav.incm.fitNonnegativePowers(obs.Rx,A,original.noisePower,settings.fitRegularization);
kkt=uav.incm.nonnegativeFitKkt(obs.Rx,A,original.noisePower,powers,settings.fitRegularization);
if kkt>settings.kktTolerance
    error('uav:incm:InvalidLocalCovfitSolution','局部CovFit功率解未满足非负KKT条件。');
end
Ai=A(:,D+1:end);pi=powers(D+1:end);
RinHat=(Ai.*pi.')*Ai'+original.noisePower*eye(cfg.derived.M);
RinHat=(RinHat+RinHat')/2;
diagnostics=struct('method',"incm_covfit_refined",'estimator',"local-covfit", ...
    'candidateManifold',A,'desiredCandidateNodes',nodes(1:D,:), ...
    'interfererCandidateNodes',nodes(D+1:end,:),'fittedPowers',powers, ...
    'desiredFittedPowers',powers(1:D),'interfererFittedPowers',pi, ...
    'noisePower',original.noisePower,'noiseSource',original.noiseSource, ...
    'refinement',refinement,'fitDiagnostics',fit,'relativeKktResidual',kkt, ...
    'usedTruth',false);
diagnostics.desiredExclusionHalfWidthDeg=uav.spectrum.desiredExclusionHalfWidth(cfg);
diagnostics.candidateGridSize=cfg.spectrum.gridSize;
diagnostics.candidateGridStep=cfg.derived.spectrumGrid.step;
diagnostics.refinementApplied=logical(applyRefinement);
diagnostics.initialCandidateNodes=initialNodes;
diagnostics.initialCandidateManifold=initialManifold;
diagnostics.peakSelection=candidateSelection;
diagnostics.loading=original.loading;
diagnostics.noiseDiagnostics=original.noiseDiagnostics;
diagnostics.effectiveSnapshotCount=original.effectiveSnapshotCount;
diagnostics.desiredCandidateFallback=original.desiredCandidateFallback;
diagnostics.commonParameters=rmfield(settings,{'radiusBins','solver'});
if ~applyRefinement
    diagnostics.method="incm_covfit";diagnostics.estimator="upa-covfit";
end
end
