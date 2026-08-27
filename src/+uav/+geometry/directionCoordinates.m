function [xi, azElRad] = directionCoordinates(uArray)
%DIRECTIONCOORDINATES Return direction cosines and azimuth/elevation in radians.

if size(uArray, 1) ~= 3 || any(~isfinite(uArray), 'all')
    error('uav:geometry:InvalidDirection', 'uArray must be finite with three rows.');
end
if any(abs(vecnorm(uArray) - 1) > 1e-10)
    error('uav:geometry:InvalidDirectionNorm', ...
        'Every array direction must have unit norm within tolerance.');
end

xi = uArray(1:2, :);
if any(sum(abs(xi).^2, 1) > 1 + 1e-10)
    error('uav:geometry:OutsideVisibleDisk', ...
        'Direction cosines must lie within the unit disk.');
end
azimuthRad = atan2(uArray(2, :), uArray(1, :));
elevationRad = atan2(uArray(3, :), hypot(uArray(1, :), uArray(2, :)));
azElRad = [azimuthRad; elevationRad];
end
