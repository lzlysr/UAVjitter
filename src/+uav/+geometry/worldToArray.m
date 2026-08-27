function uArray = worldToArray(R, uWorld)
%WORLDTOARRAY Express world-frame directions in the instantaneous array frame.

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

uArray = R.' * uWorld;
end
