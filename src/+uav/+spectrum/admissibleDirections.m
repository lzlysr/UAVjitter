function mask = admissibleDirections(nodes, cfg)
%ADMISSIBLEDIRECTIONS Test continuous lower-hemisphere INCM support nodes.

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
z = -sqrt(max(0, 1 - radialSquared));
directions = [real(nodes), z];
nominalDirection = cfg.derived.nominalArrayUnitVectors(:, 1);
cosine = directions * nominalDirection;
angularDistance = acos(max(-1, min(1, real(cosine))));
outsideSoi = angularDistance > ...
    deg2rad(cfg.spectrum.desiredExclusionHalfWidthDeg);
mask = visible & outsideSoi;
end
