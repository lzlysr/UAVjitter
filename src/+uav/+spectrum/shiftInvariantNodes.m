function [nodes, powers] = shiftInvariantNodes(Rx, noisePower, cfg)
%SHIFTINVARIANTNODES 从当前SCM的平移不变子空间提取并细化二维候选。
%   用平移不变方程初始化联合子空间拟合。候选的跨周期归属由调用方决定。

nx = cfg.array.Nx;
ny = cfg.array.Ny;
M = nx*ny;
settings = cfg.incm.temporalCovfit;
index = reshape(1:M, ny, nx);
xFirst = index(:, 1:end-1);
xSecond = index(:, 2:end);
yFirst = index(1:end-1, :);
ySecond = index(2:end, :);

[vectors, eigenvalues] = eig((Rx + Rx')/2, 'vector');
[eigenvalues, order] = sort(real(eigenvalues), 'descend');
maximumRank = min([numel(xFirst), numel(yFirst), cfg.incm.localCovfit.maximumCandidateCount]);
signalRank = min(maximumRank, nnz(eigenvalues > settings.minimumModeToNoise * noisePower));
nodes = zeros(0, 2);
powers = zeros(0, 1);
if signalRank == 0
    return;
end

U = vectors(:, order(1:signalRank));
xShift = U(xFirst(:), :) \ U(xSecond(:), :);
yShift = U(yFirst(:), :) \ U(ySecond(:), :);
% 固定复系数将两轴的共同特征基分离，也允许某一轴出现相同空间频率。
[basis, ~] = eig(xShift + exp(1i*pi/4) * yShift);
if rcond(basis) <= sqrt(eps)
    return;
end
phase = angle([diag(basis \ (xShift*basis)), diag(basis \ (yShift*basis))]);
nodes = phase ./ (2*pi*cfg.array.spacingOverLambda);
eigenvalues = eigenvalues(1:signalRank).';
% 强且可靠的子空间方向对联合拟合有较大权重。
weightedBasis = U .* ((eigenvalues-noisePower) ./ sqrt(eigenvalues));
phaseSlope = 2*pi*[kron((0:nx-1)',ones(ny,1)), repmat((0:ny-1)',nx,1)] ...
    .* cfg.array.spacingOverLambda;
nodes = refineNodes(nodes, phaseSlope, weightedBasis, settings.modeRefinementIterations);
% 边缘处的轻微估计误差可投回圆盘，明显越界的模态直接拒绝。
radius = vecnorm(nodes, 2, 2);
keep = radius <= 1 + settings.mergeDistance;
nodes = nodes(keep, :) ./ max(1, radius(keep));
if isempty(nodes)
    return;
end
A = exp(1i*phaseSlope*nodes.');
powers = uav.spectrum.candidateAssociationPowers(A,Rx,noisePower);
if isempty(powers), nodes = zeros(0,2); end
end

function nodes = refineNodes(nodes, phaseSlope, Y, iterationCount)
% 直接减小当前有效子空间中无法由候选流形解释的能量。
for iteration = 1:iterationCount
    A = exp(1i*phaseSlope*nodes.');
    if rcond(A'*A) <= sqrt(eps)
        return;
    end
    coefficients = A \ Y;
    residual = Y-A*coefficients;
    derivative = zeros(numel(residual),numel(nodes),'like',Y);
    for source = 1:size(nodes,1)
        tangent = 1i*phaseSlope .* A(:,source);
        tangent = tangent - A*(A \ tangent);
        for axisIndex = 1:2
            block = tangent(:,axisIndex)*coefficients(source,:);
            derivative(:,source+(axisIndex-1)*size(nodes,1)) = block(:);
        end
    end
    jacobian = [real(derivative);imag(derivative)];
    if rcond(jacobian'*jacobian) <= sqrt(eps)
        return;
    end
    step = reshape(jacobian \ [real(residual(:));imag(residual(:))],size(nodes));
    candidate = nodes+step;
    candidateA = exp(1i*phaseSlope*candidate.');
    if rcond(candidateA'*candidateA) <= sqrt(eps) || ...
            norm(Y-candidateA*(candidateA\Y),'fro') >= norm(residual,'fro')
        return;
    end
    nodes = candidate;
end
end
