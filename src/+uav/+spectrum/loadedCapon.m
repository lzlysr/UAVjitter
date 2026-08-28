function [powerMap, diagnostics] = loadedCapon(Rx, cfg)
%LOADEDCAPON Evaluate loaded Capon power over a visible direction-cosine grid.

validateInputs(Rx, cfg);
M = cfg.derived.M;
Rx = (Rx + Rx') / 2;
loading = cfg.spectrum.diagonalLoading * covarianceScale(Rx);
loadedCovariance = Rx + loading * eye(M);
grid = cfg.derived.spectrumGrid;
manifold = grid.manifold;
response = loadedCovariance \ manifold;
if any(~isfinite(response), 'all')
    response = pinv(loadedCovariance) * manifold;
end
denominator = real(sum(conj(manifold) .* response, 1));
tolerance = eps(max(1, sqrt(M) * vecnorm(response, 2, 1)));
if any(denominator <= tolerance)
    error('uav:spectrum:DegenerateCovariance', ...
        'Loaded covariance has a nonpositive Capon denominator.');
end
powerMap = nan(size(grid.xiX));
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
