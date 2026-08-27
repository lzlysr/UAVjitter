function R = rotationZYX(etaRad)
%ROTATIONZYX Return the roll-pitch-yaw rotation from array to world axes.

etaRad = etaRad(:);
if numel(etaRad) ~= 3 || any(~isfinite(etaRad))
    error('uav:geometry:InvalidAttitude', 'etaRad must contain three finite radians.');
end

roll = etaRad(1);
pitch = etaRad(2);
yaw = etaRad(3);
Rx = [1, 0, 0; 0, cos(roll), -sin(roll); 0, sin(roll), cos(roll)];
Ry = [cos(pitch), 0, sin(pitch); 0, 1, 0; -sin(pitch), 0, cos(pitch)];
Rz = [cos(yaw), -sin(yaw), 0; sin(yaw), cos(yaw), 0; 0, 0, 1];
R = Rz * Ry * Rx;
end
