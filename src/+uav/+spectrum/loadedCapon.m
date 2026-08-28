function [powerMap, diagnostics] = loadedCapon(Rx, cfg)
%LOADEDCAPON 在二维方向余弦可见域计算加载Capon空间谱。
%   Rx       - M×M当前周期样本协方差；
%   powerMap - gridSize×gridSize功率图，可见圆盘外为NaN；
%   diagnostics包含网格、加载量和实际求解使用的协方差。
%
%   Capon谱为P(xi)=1/(a(xi)^H Rload^{-1} a(xi))。网格流形在配置
%   校验阶段预计算，此处把所有可见节点组成M×N矩阵并一次批量求解。

validateInputs(Rx, cfg);
M = cfg.derived.M;
Rx = (Rx + Rx') / 2;
% 加载量随平均阵元功率缩放，避免信号整体缩放时参数失去含义。
loading = cfg.spectrum.diagonalLoading * covarianceScale(Rx);
loadedCovariance = Rx + loading * eye(M);
grid = cfg.derived.spectrumGrid;
manifold = grid.manifold;
response = loadedCovariance \ manifold;
if any(~isfinite(response), 'all')
    response = pinv(loadedCovariance) * manifold;
end
% 每一列实现a^H*(Rload\a)，避免显式求逆矩阵。
denominator = real(sum(conj(manifold) .* response, 1));
tolerance = eps(max(1, sqrt(M) * vecnorm(response, 2, 1)));
if any(denominator <= tolerance)
    error('uav:spectrum:DegenerateCovariance', ...
        'Loaded covariance has a nonpositive Capon denominator.');
end
powerMap = nan(size(grid.xiX));
% 使用与缓存流形相同的线性索引顺序写回二维功率图。
powerMap(grid.visibleLinearIndices) = 1 ./ denominator;
diagnostics = struct( ...
    'xiX', grid.xiX, ...
    'xiY', grid.xiY, ...
    'visibleMask', grid.visibleMask, ...
    'visibleLinearIndices', grid.visibleLinearIndices, ...
    'cellArea', grid.cellArea, ...
    'usedCachedManifold', true, ...
    'loading', loading, ...
    'loadedCovariance', loadedCovariance);
end

function validateInputs(Rx, cfg)
if ~isstruct(cfg) || ~isfield(cfg, 'derived') || ~isfield(cfg.derived, 'spectrumGrid') || ...
        ~isfield(cfg, 'array') || ~isfield(cfg, 'spectrum')
    error('uav:spectrum:InvalidConfiguration', 'cfg is missing spectrum fields.');
end
M = cfg.derived.M;
if ~isequal(size(Rx), [M, M]) || any(~isfinite(Rx), 'all')
    error('uav:spectrum:InvalidCovariance', ...
        'Rx must be finite and match the configured array size.');
end
end

function scale = covarianceScale(Rx)
scale = max(real(trace(Rx)) / size(Rx, 1), eps);
end
