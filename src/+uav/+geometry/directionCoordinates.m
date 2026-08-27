function [xi, azElRad] = directionCoordinates(uArray)
%DIRECTIONCOORDINATES Return direction cosines and azimuth/elevation in radians.

if size(uArray, 1) ~= 3 || any(~isfinite(uArray), 'all')
    error('uav:geometry:InvalidDirection', 'uArray must be finite with three rows.');
end

xi = uArray(1:2, :);
azimuthRad = atan2(uArray(2, :), uArray(1, :));
elevationRad = atan2(uArray(3, :), hypot(uArray(1, :), uArray(2, :)));
azElRad = [azimuthRad; elevationRad];
end
