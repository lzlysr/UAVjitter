function [nodes,powers] = fullCaponNodes(obs,noisePower,cfg)
%FULLCAPONNODES 对可见域内的Capon局部峰做连续细化。
%   峰在区内还是区外均可先成为候选，由调用方结合历史决定干扰身份。
%   不使用历史来约束本函数的峰值位置。
original = struct('desiredCandidateNodes',zeros(0,2));
[~,~,nodes] = uav.spectrum.localPeakCandidates(obs.Rx,original,cfg,true);
[A,nodes] = uav.spectrum.refineInterfererCandidates(obs,nodes,0,cfg,true);
% 不同网格峰可能收敛至同一位置，保留原谱值排序更靠前的一项。
keep = false(size(nodes,1),1);
for index = 1:size(nodes,1)
    keep(index) = all(vecnorm(nodes(keep,:)-nodes(index,:),2,2) > ...
        cfg.incm.temporalCovfit.mergeDistance);
end
nodes = nodes(keep,:);
powers = uav.spectrum.candidateAssociationPowers(A(:,keep),obs.Rx,noisePower);
if isempty(powers), nodes = zeros(0,2); end
end
