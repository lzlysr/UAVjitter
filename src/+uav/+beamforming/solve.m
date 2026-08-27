function [w, diagnostics] = solve(method, R, aRef, cfg)
%SOLVE Compute MVDR or diagonally loaded MVDR weights.

method = lower(string(method));
validateInputs(R, aRef);
R = (R + R') / 2;
M = size(R, 1);

switch method
    case "mvdr"
        loading = minimumNumericalLoading(R, cfg);
    case "dlmvdr"
        loading = cfg.beamforming.diagonalLoading * covarianceScale(R);
    otherwise
        error('uav:beamforming:UnsupportedMethod', ...
            'Task 3 supports only "mvdr" and "dlmvdr".');
end

loadedCovariance = R + loading * eye(M);
rawWeight = loadedCovariance \ aRef;
if any(~isfinite(rawWeight)) || norm(rawWeight) == 0
    rawWeight = pinv(loadedCovariance) * aRef;
end
normalizer = aRef' * rawWeight;
if ~isfinite(normalizer) || abs(normalizer) <= eps(max(1, norm(aRef) * norm(rawWeight)))
    error('uav:beamforming:DegenerateReference', ...
        'The reference steering vector cannot be normalized against R.');
end
w = rawWeight / normalizer;
diagnostics = struct( ...
    'method', method, ...
    'loading', loading, ...
    'loadedCovariance', loadedCovariance, ...
    'reciprocalCondition', rcond(loadedCovariance), ...
    'normalizer', normalizer);
end

function validateInputs(R, aRef)
if ~ismatrix(R) || size(R, 1) ~= size(R, 2) || any(~isfinite(R), 'all')
    error('uav:beamforming:InvalidCovariance', ...
        'R must be a finite square covariance matrix.');
end
if ~isvector(aRef) || numel(aRef) ~= size(R, 1) || any(~isfinite(aRef))
    error('uav:beamforming:InvalidReference', ...
        'aRef must be a finite vector matching R.');
end
end

function scale = covarianceScale(R)
scale = max(real(trace(R)) / size(R, 1), eps);
end

function loading = minimumNumericalLoading(R, cfg)
minimumReciprocalCondition = cfg.beamforming.minimumReciprocalCondition;
if rcond(R) >= minimumReciprocalCondition
    loading = 0;
    return;
end

eigenvalues = real(eig(R));
scale = covarianceScale(R);
largestEigenvalue = max(max(eigenvalues), scale);
smallestEigenvalue = min(eigenvalues);
eigenvalueLoading = max(0, (minimumReciprocalCondition * largestEigenvalue - ...
    smallestEigenvalue) / (1 - minimumReciprocalCondition));
oneNormTarget = minimumReciprocalCondition * sqrt(size(R, 1));
oneNormLoading = oneNormTarget * max(norm(R, 1), scale) / (1 - oneNormTarget);
loading = max(eigenvalueLoading, oneNormLoading);
end
