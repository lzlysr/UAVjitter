function [A, desiredCount, nodes, diagnostics] = localPeakCandidates(Rx, original, cfg)
%LOCALPEAKCANDIDATES 只改变干扰候选：二维局部极大值，而不是最强相邻网格单元。
%   期望字典沿用原CovFit；最终波束参考始终由实验层固定为名义导向。
%   不使用真实干扰数或方向；阈值、上限、加载沿用原配置。
local = cfg;
local.spectrum.diagonalLoading=cfg.incm.staticCovfit.loadingRatio;
local.spectrum.peakThresholdRatio=cfg.incm.staticCovfit.peakThresholdRatio;
local.spectrum.peakNmsRadiusBins=cfg.incm.staticCovfit.peakNmsRadiusBins;
[power,s] = uav.spectrum.loadedCapon(Rx,local);
map=power; map(~s.visibleMask)=-Inf;
isPeak=s.visibleMask;
% 八邻域局部极大值检测；不跨越网格边界环绕比较。
for dr=-1:1
    for dc=-1:1
        if dr==0 && dc==0,continue;end
        neighbor=-inf(size(map));
        rows=max(1,1+dr):min(size(map,1),size(map,1)+dr);
        cols=max(1,1+dc):min(size(map,2),size(map,2)+dc);
        neighbor(rows,cols)=map(rows-dr,cols-dc);
        isPeak=isPeak & map>=neighbor;
    end
end
allNodes=[s.xiX(:),s.xiY(:)];
outside=reshape(uav.spectrum.admissibleDirections(allNodes,cfg),size(map));
[interferers,selection] = uav.spectrum.selectPeaks(power,s,isPeak & outside,local);
cap=min(cfg.incm.staticCovfit.maximumCandidateCount,cfg.derived.M-1);
interferers=interferers(1:min(size(interferers,1),cap),:);
desired=original.desiredCandidateNodes;
nodes=[desired;interferers]; desiredCount=size(desired,1);
A=uav.array.steeringUPA(nodes.',cfg.array.Nx,cfg.array.Ny,cfg.array.spacingOverLambda);
diagnostics=struct('selection',selection,'interfererCount',size(interferers,1),'usedTruth',false);
end
