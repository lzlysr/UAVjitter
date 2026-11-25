function [RinHat, diagnostics] = subspaceCovfit(obs,cfg)
%SUBSPACECOVFIT 仅用二维平移不变子空间定向，并联合拟合点源功率。
%   单周期用固定期望角区划分候选。保留全域候选供时序版做历史关联。
if min(cfg.array.Nx,cfg.array.Ny) < 2
    error('uav:incm:SubspaceCovfitRequiresUpa','二维纯子空间CovFit要求UPA两轴均有相邻阵元。');
end
L = 0;
if isfield(obs,'X'), L = size(obs.X,2); end
[noise,noiseSource] = uav.incm.covfitNoise(obs.Rx,L,cfg);
[nodes,associationPowers] = uav.spectrum.shiftInvariantNodes(obs.Rx,noise,cfg);
outside = uav.spectrum.admissibleDirections(nodes,cfg);
nominal = cfg.derived.nominalDirectionCosines(:,1).';
desired = nominal;
inside = nodes(~outside,:);
% 以三维夹角最近的区内候选解释期望分量，避免默认选择区内最强来波。
if ~isempty(inside)
    directions = [inside, -sqrt(max(0,1-sum(inside.^2,2)))];
    [~,nearest] = max(directions*cfg.derived.nominalArrayUnitVectors(:,1));
    node = inside(nearest,:);
    if norm(node-nominal,inf) > cfg.spectrum.quadratureDedupTolerance
        desired = [desired;node];
    end
end
interferers = nodes(outside,:);
D = size(desired,1);
A = uav.array.steeringUPA([desired;interferers].',cfg.array.Nx,cfg.array.Ny, ...
    cfg.array.spacingOverLambda);
settings = cfg.incm.localCovfit;
[p,fit] = uav.incm.fitNonnegativePowers(obs.Rx,A,noise,settings.fitRegularization);
kkt = uav.incm.nonnegativeFitKkt(obs.Rx,A,noise,p,settings.fitRegularization);
if kkt > settings.kktTolerance
    error('uav:incm:InvalidSubspaceCovfitSolution','纯子空间CovFit功率解未满足非负KKT条件。');
end
Ai = A(:,D+1:end);
interfererPowers = p(D+1:end);
RinHat = (Ai.*interfererPowers.')*Ai' + noise*eye(cfg.derived.M);
RinHat = (RinHat+RinHat')/2;
diagnostics = struct('method',"incm_covfit_refined_temporal_pure_subspace",'usedTruth',false, ...
    'temporalApplied',false,'noisePower',noise,'noiseSource',noiseSource, ...
    'candidateManifold',A,'desiredCandidateNodes',desired, ...
    'interfererCandidateNodes',interferers,'proposedCandidateNodes',nodes, ...
    'proposedAssociationPowers',associationPowers,'fittedPowers',p, ...
    'desiredFittedPowers',p(1:D),'interfererFittedPowers',interfererPowers, ...
    'fitDiagnostics',fit,'relativeKktResidual',kkt);
end
