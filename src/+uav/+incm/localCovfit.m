function [RinHat,diagnostics] = localCovfit(obs,cfg)
%LOCALCOVFIT Capon局部干扰峰细化与非负协方差拟合，固定最终期望参考。
%   移植自研究分支1477a46的refined_nnls0，不是Le全分区版本。
%   流程：低加载谱→局部峰→有界连续干扰方向细化→联合NNLS→剔除期望组。
%   只细化干扰节点，期望拟合字典不变；不读取任何真实方向或真实INCM。
settings=cfg.incm.localCovfit;
front=cfg;
% 复用旧接口获取期望候选、噪声及网格信息。其额外旧拟合结果不参与新INCM。
front.incm.staticCovfit.loadingRatio=settings.loadingRatio;
front.incm.staticCovfit.peakThresholdRatio=settings.peakThresholdRatio;
front.incm.staticCovfit.peakNmsRadiusBins=settings.peakNmsRadiusBins;
front.incm.staticCovfit.maximumCandidateCount=settings.maximumCandidateCount;
[~,original]=uav.incm.upaCovfit(obs,front);
[A,nodes,refinement]=uav.spectrum.refineInterfererCandidates(obs,original,front);
[powers,fit]=uav.incm.fitNonnegativePowers(obs.Rx,A,original.noisePower,settings.fitRegularization);
kkt=uav.incm.nonnegativeFitKkt(obs.Rx,A,original.noisePower,powers,settings.fitRegularization);
if kkt>settings.kktTolerance
    error('uav:incm:InvalidLocalCovfitSolution','局部CovFit功率解未满足非负KKT条件。');
end
D=refinement.desiredCount;Ai=A(:,D+1:end);pi=powers(D+1:end);
RinHat=(Ai.*pi.')*Ai'+original.noisePower*eye(cfg.derived.M);
RinHat=(RinHat+RinHat')/2;
diagnostics=struct('method',"CovFit_local",'estimator',"local-covfit", ...
    'candidateManifold',A,'desiredCandidateNodes',nodes(1:D,:), ...
    'interfererCandidateNodes',nodes(D+1:end,:),'fittedPowers',powers, ...
    'desiredFittedPowers',powers(1:D),'interfererFittedPowers',pi, ...
    'noisePower',original.noisePower,'noiseSource',original.noiseSource, ...
    'refinement',refinement,'fitDiagnostics',fit,'relativeKktResidual',kkt, ...
    'usedTruth',false);
end
