function uArray = worldToArray(R, uWorld)
%WORLDTOARRAY Express world-frame directions in the instantaneous array frame.

if ~isequal(size(R), [3, 3]) || any(~isfinite(R), 'all')
    error('uav:geometry:InvalidRotation', 'R must be a finite 3-by-3 rotation.');
end
if size(uWorld, 1) ~= 3 || any(~isfinite(uWorld), 'all')
    error('uav:geometry:InvalidDirection', 'uWorld must be finite with three rows.');
end

uArray = R.' * uWorld;
end
