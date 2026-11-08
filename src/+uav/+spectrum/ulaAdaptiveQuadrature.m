function quadrature = ulaAdaptiveQuadrature(power, spectrum, mask, peaks, cfg)
%ULAADAPTIVEQUADRATURE 在ULA等角网格上进行一维局部细分及粗块求积。
%   节点仍输出[xi_x,0]，但积分权重始终是弧度长度dtheta。
%   原始允许单元互斥分组；细单元权重和、粗块权重和均守恒。
theta = spectrum.thetaRad;
step = spectrum.step;
fineMask = false(size(mask));
radius = cfg.spectrum.quadratureFineRadiusBins;
for p = reshape(peaks,1,[])
    fineMask(max(1,p-radius):min(numel(mask),p+radius)) = true;
end
fineMask = fineMask & mask;
nodes = zeros(0,2); weights = zeros(0,1); values = zeros(0,1);
refinement = cfg.spectrum.quadratureRefinement;
offsets = (((1:refinement)-0.5)/refinement-0.5)*step;
for index = reshape(find(fineMask),1,[])
    angles = theta(index)+offsets;
    candidates = [sin(angles(:)),zeros(numel(angles),1)];
    keep = uav.spectrum.admissibleDirections(candidates,cfg);
    candidates = candidates(keep,:); angles = angles(keep);
    if isempty(angles)
        angles = theta(index); candidates = [sin(angles),0];
    end
    interpolated = interp1(theta,power,angles,'linear');
    interpolated(~isfinite(interpolated)) = power(index);
    nodes = [nodes;candidates]; %#ok<AGROW>
    weights = [weights;step/numel(angles)*ones(numel(angles),1)]; %#ok<AGROW>
    values = [values;max(real(interpolated(:)),0)]; %#ok<AGROW>
end
fineCount = size(nodes,1);
indices = find(mask & ~fineMask);
groups = floor((indices-1)/cfg.spectrum.quadratureCoarseStride);
for group = reshape(unique(groups),1,[])
    members = indices(groups==group);
    [~,nearest] = min(abs(theta(members)-mean(theta(members))));
    nodes(end+1,:) = [sin(theta(members(nearest))),0]; %#ok<AGROW>
    weights(end+1,1) = numel(members)*step; %#ok<AGROW>
    values(end+1,1) = mean(power(members)); %#ok<AGROW>
end
coarseCount = size(nodes,1)-fineCount;
if isempty(nodes)
    merged=struct('deduplicatedNodeCount',0);
else
    [nodes,weights,values,merged] = uav.spectrum.mergeQuadratureNodes( ...
        nodes,weights,values,cfg.spectrum.quadratureDedupTolerance);
end
quadrature = struct('nodes',nodes,'weights',weights,'values',values, ...
    'diagnostics',struct('mode',"adaptive-ULA",'measure',"dtheta", ...
    'fineCellCount',nnz(fineMask),'coarseCellCount',numel(indices), ...
    'fineNodeCount',fineCount,'coarseNodeCount',coarseCount, ...
    'deduplicatedNodeCount',merged.deduplicatedNodeCount, ...
    'representedArea',sum(weights),'fineMask',fineMask));
end
