function [powerMap, diagnostics] = loadedCapon(Rx, cfg)
%LOADEDCAPON Evaluate loaded Capon power over a visible direction-cosine grid.

validateInputs(Rx, cfg);
M = cfg.derived.M;
Rx = (Rx + Rx') / 2;
loading = cfg.spectrum.diagonalLoading * covarianceScale(Rx);
loadedCovariance = Rx + loading * eye(M);
gridAxis = linspace(-1, 1, cfg.spectrum.gridSize);
[xiX, xiY] = meshgrid(gridAxis, gridAxis);
visibleMask = xiX.^2 + xiY.^2 <= 1;
powerMap = nan(size(xiX));

for linearIndex = find(visibleMask).'
    a = uav.array.steeringUPA([xiX(linearIndex); xiY(linearIndex)], ...
        cfg.array.Nx, cfg.array.Ny, cfg.array.spacingOverLambda);
    response = loadedCovariance \ a;
    if any(~isfinite(response))
        response = pinv(loadedCovariance) * a;
    end
    denominator = real(a' * response);
    if denominator <= eps(max(1, norm(a) * norm(response)))
        error('uav:spectrum:DegenerateCovariance', ...
            'Loaded covariance has a nonpositive Capon denominator.');
    end
    powerMap(linearIndex) = 1 / denominator;
end

if cfg.spectrum.gridSize == 1
    cellArea = 4;
else
    gridStep = gridAxis(2) - gridAxis(1);
    cellArea = gridStep^2;
end
diagnostics = struct( ...
    'xiX', xiX, ...
    'xiY', xiY, ...
    'visibleMask', visibleMask, ...
    'cellArea', cellArea, ...
    'loading', loading, ...
    'loadedCovariance', loadedCovariance);
end

function validateInputs(Rx, cfg)
if ~isstruct(cfg) || ~isfield(cfg, 'derived') || ~isfield(cfg, 'array') || ...
        ~isfield(cfg, 'spectrum')
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
