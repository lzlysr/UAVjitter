function [refinedNodes, diagnostics] = refineCaponNodesLe(nodes, covariance, insideDesired, cfg)
%   移植自Gu–Leshem研究分支d666399；仅适配配置路径及函数名。
%REFINECAPONNODES 在原先验分区内局部细化Capon候选，减小离网误差。
%   本函数只用传入协方差、候选起点、阵列与固定期望角域，不读取真值。
%   最小化a(xi)^H*R^{-1}*a(xi)等价于最大化Capon谱。
%   该方向目标一般非凸，fmincon仅给局部解，不保证发现所有真实源。
%   输入R必须已按调用方约定处理，本函数不额外加加载或使用伪逆。
if size(nodes, 2) ~= 2 || any(~isfinite(nodes), 'all') || ...
        ~isequal(size(covariance), [cfg.derived.M, cfg.derived.M]) || ...
        any(~isfinite(covariance), 'all') || rcond(covariance) <= eps
    error('uav:spectrum:InvalidRefinementInput', '候选或细化协方差尺寸非法、非有限或奇异。');
end
settings = cfg.experiments.linYueStudy.leRefinement;
options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'off', ...
    'SpecifyObjectiveGradient', true, 'MaxIterations', settings.maxIterations, ...
    'OptimalityTolerance', settings.optimalityTolerance, ...
    'ConstraintTolerance', settings.constraintTolerance, ...
    'StepTolerance', settings.stepTolerance);
% 阵元顺序必须与vx kron vy一致：y索引先变化。
indexX = repelem((0:cfg.array.Nx-1).', cfg.array.Ny);
indexY = repmat((0:cfg.array.Ny-1).', cfg.array.Nx, 1);
phaseSlopes = 2*pi * [cfg.array.spacingOverLambda(1)*indexX, ...
    cfg.array.spacingOverLambda(2)*indexY];
covariance = (covariance + covariance') / 2;
refinedNodes = nodes;
exitflags = zeros(size(nodes, 1), 1);
accepted = false(size(nodes, 1), 1);
for index = 1:size(nodes, 1)
    initial = nodes(index, :).';
    [solution, ~, exitflag] = fmincon(@objective, initial, ...
        [], [], [], [], -ones(2, 1), ones(2, 1), @constraints, options);
    % 求解容差可能使节点略超圆盘；投回物理边界后重新检查目标与角域。
    if norm(solution) > 1, solution = solution / norm(solution); end
    exitflags(index) = exitflag;
    violations = constraints(solution);
    if exitflag > 0 && all(isfinite(solution)) && ...
            all(violations <= settings.constraintTolerance) && ...
            objective(solution) <= objective(initial) + settings.optimalityTolerance
        refinedNodes(index, :) = solution.';
        accepted(index) = true;
    end
end
diagnostics = struct('exitflags', exitflags, 'accepted', accepted, ...
    'initialNodes', nodes, 'refinedNodes', refinedNodes, 'usedTruth', false);

    function [value, gradient] = objective(xi)
        % 直接构造同一UPA流形，给出二维实参数下的解析梯度。
        a = exp(1i * phaseSlopes * xi);
        solved = covariance \ a;
        value = real(a' * solved);
        derivative = 1i * phaseSlopes .* a;
        gradient = 2 * real(derivative' * solved);
    end
    function [inequality, equality] = constraints(xi)
        % 可见圆盘和固定三维角域约束，不随当前真实姿态更新。
        radiusSquared = sum(xi.^2);
        direction = [xi; -sqrt(max(0, 1-radiusSquared))];
        cosine = cfg.derived.nominalArrayUnitVectors(:, 1)' * direction;
        boundary = cosd(cfg.spectrum.desiredExclusionHalfWidthDeg);
        if insideDesired
            regionConstraint = boundary - cosine;
        else
            regionConstraint = cosine - boundary;
        end
        inequality = [radiusSquared-1; regionConstraint];
        equality = [];
    end
end
