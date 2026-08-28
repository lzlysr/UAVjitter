function mask = admissibleDirections(nodes, cfg)
%ADMISSIBLEDIRECTIONS 判断连续方向节点能否参与干扰谱积分。
%   nodes为N×2方向余弦[xi_x,xi_y]。节点必须同时满足：
%   1) 位于单位可见圆盘；2) 位于名义期望方向排除角域之外。
%   使用连续三维夹角而不是矩形网格距离，保证粗/细求积节点规则一致。

if ~isnumeric(nodes) || size(nodes, 2) ~= 2 || any(~isfinite(nodes), 'all')
    error('uav:spectrum:InvalidDirectionNodes', ...
        'nodes must be a finite N-by-2 direction-cosine matrix.');
end
if ~isstruct(cfg) || ~isfield(cfg, 'spectrum') || ~isfield(cfg, 'derived') || ...
        ~isfield(cfg.spectrum, 'desiredExclusionHalfWidthDeg') || ...
        ~isfield(cfg.derived, 'nominalArrayUnitVectors')
    error('uav:spectrum:InvalidConfiguration', ...
        'cfg is missing continuous admissibility fields.');
end

radialSquared = sum(real(nodes).^2, 2);
visible = radialSquared <= 1;
% 当前场景信号来自阵列下半球，取负z分量恢复三维单位方向。
z = -sqrt(max(0, 1 - radialSquared));
directions = [real(nodes), z];
nominalDirection = cfg.derived.nominalArrayUnitVectors(:, 1);
cosine = directions * nominalDirection;
angularDistance = acos(max(-1, min(1, real(cosine))));
outsideSoi = angularDistance > ...
    deg2rad(cfg.spectrum.desiredExclusionHalfWidthDeg);
mask = visible & outsideSoi;
end
