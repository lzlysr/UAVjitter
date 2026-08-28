function uArray = worldToArray(R, uWorld)
%WORLDTOARRAY 将世界坐标系方向转换到当前阵列坐标系。
%   R      - 3×3正旋转矩阵，满足v_world=R*v_array；
%   uWorld - 3×Q世界坐标单位方向，每列对应一个信号源；
%   uArray - 3×Q阵列坐标单位方向。
%
%   因为R是正交矩阵，其逆矩阵等于转置，所以反向坐标变换使用R'。

if ~isequal(size(R), [3, 3]) || any(~isfinite(R), 'all') || ...
        norm(R' * R - eye(3), 'fro') > 1e-10 || abs(det(R) - 1) > 1e-10
    error('uav:geometry:InvalidRotation', ...
        'R must be a finite proper rotation within tolerance.');
end
if size(uWorld, 1) ~= 3 || any(~isfinite(uWorld), 'all')
    error('uav:geometry:InvalidDirection', 'uWorld must be finite with three rows.');
end
if any(abs(vecnorm(uWorld) - 1) > 1e-10)
    error('uav:geometry:InvalidDirectionNorm', ...
        'Every world direction must have unit norm within tolerance.');
end

% 只改变坐标表示，不改变每列向量的长度。
uArray = R.' * uWorld;
end
