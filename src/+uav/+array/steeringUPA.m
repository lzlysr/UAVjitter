function a = steeringUPA(xi, Nx, Ny, dOverLambda)
%STEERINGUPA 生成矩形均匀平面阵列的窄带导向矢量。
%   xi          - 2×Q方向余弦，每列为[xi_x;xi_y]；
%   Nx, Ny      - 阵列x/y方向的阵元数；
%   dOverLambda - [dx/lambda,dy/lambda]，也可输入相同间距的标量；
%   a           - (Nx*Ny)×Q复导向矢量矩阵。
%
%   阵元顺序固定为kron(vx,vy)，即y方向索引变化最快。每个阵元幅度为1，
%   所以每列导向矢量的2范数为sqrt(Nx*Ny)。本函数支持一次输入多个方向。

if ~isnumeric(xi) || size(xi, 1) ~= 2 || any(~isfinite(xi), 'all')
    error('uav:array:InvalidDirectionCosines', 'xi must be finite with two rows.');
end
if any(sum(abs(xi).^2, 1) > 1 + 1e-10)
    error('uav:array:OutsideVisibleDisk', ...
        'Each direction-cosine column must lie within the unit disk.');
end
if ~isscalar(Nx) || Nx <= 0 || Nx ~= floor(Nx) || ...
        ~isscalar(Ny) || Ny <= 0 || Ny ~= floor(Ny)
    error('uav:array:InvalidDimensions', 'Nx and Ny must be positive integers.');
end
if isscalar(dOverLambda)
    dOverLambda = [dOverLambda, dOverLambda];
end
if ~isnumeric(dOverLambda) || numel(dOverLambda) ~= 2 || ...
        any(~isfinite(dOverLambda)) || any(dOverLambda <= 0)
    error('uav:array:InvalidSpacing', 'dOverLambda must contain two positive finite spacings.');
end
dOverLambda = reshape(dOverLambda, 1, []);

% x和y两个一维均匀线阵的相位响应。
vx = exp(1i * 2 * pi * dOverLambda(1) * xi(1, :) .* (0:Nx - 1).');
vy = exp(1i * 2 * pi * dOverLambda(2) * xi(2, :) .* (0:Ny - 1).');
% 通过Kronecker积组合成矩形UPA响应，并显式保持约定的传感器顺序。
a = zeros(Nx * Ny, size(xi, 2));
for sourceIndex = 1:size(xi, 2)
    a(:, sourceIndex) = kron(vx(:, sourceIndex), vy(:, sourceIndex));
end
end
