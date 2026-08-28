function [xi, azElRad] = directionCoordinates(uArray)
%DIRECTIONCOORDINATES 从三维单位方向提取阵列建模所需坐标。
%   uArray - 3×Q阵列坐标单位方向；
%   xi     - 2×Q方向余弦[x分量;y分量]，直接进入UPA导向矢量；
%   azElRad- 2×Q[方位角;高度角]，主要用于诊断和绘图。

if size(uArray, 1) ~= 3 || any(~isfinite(uArray), 'all')
    error('uav:geometry:InvalidDirection', 'uArray must be finite with three rows.');
end
if any(abs(vecnorm(uArray) - 1) > 1e-10)
    error('uav:geometry:InvalidDirectionNorm', ...
        'Every array direction must have unit norm within tolerance.');
end

% UPA位于阵列x-y平面，因此相位只由方向的x、y分量决定。
xi = uArray(1:2, :);
if any(sum(abs(xi).^2, 1) > 1 + 1e-10)
    error('uav:geometry:OutsideVisibleDisk', ...
        'Direction cosines must lie within the unit disk.');
end
azimuthRad = atan2(uArray(2, :), uArray(1, :));
elevationRad = atan2(uArray(3, :), hypot(uArray(1, :), uArray(2, :)));
azElRad = [azimuthRad; elevationRad];
end
