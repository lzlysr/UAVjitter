function [nodes,powers] = hybridCaponNodes(obs,noisePower,cfg)
%HYBRIDCAPONNODES 期望排除区内用子空间候选，区外用Capon候选。
%   先独立细化来向，再按已有的连续三维夹角规则划分。
%   边界两侧候选重叠时优先保留区内子空间点，联合重估关联功率。
subspace = uav.spectrum.shiftInvariantNodes(obs.Rx,noisePower,cfg);
capon = uav.spectrum.fullCaponNodes(obs,noisePower,cfg);
inside = subspace(~uav.spectrum.admissibleDirections(subspace,cfg),:);
outside = capon(uav.spectrum.admissibleDirections(capon,cfg),:);

distance = hypot(outside(:,1)-inside(:,1).',outside(:,2)-inside(:,2).');
outside = outside(all(distance > cfg.incm.temporalCovfit.mergeDistance,2),:);
capacity = min(cfg.derived.M-1,cfg.incm.localCovfit.maximumCandidateCount);
inside = inside(1:min(size(inside,1),capacity),:);
remaining = min(size(outside,1),capacity-size(inside,1));
nodes = [inside;outside(1:remaining,:)];
A = uav.array.steeringUPA(nodes.',cfg.array.Nx,cfg.array.Ny,cfg.array.spacingOverLambda);
powers = uav.spectrum.candidateAssociationPowers(A,obs.Rx,noisePower);
if isempty(powers), nodes = zeros(0,2); end
end
