function [w, diagnostics] = solve(method, R, aRef, cfg)
%SOLVE Compute MVDR, diagonally loaded MVDR, or robust common-INCM weights.

method = lower(string(method));
validateInputs(R, aRef, cfg);
R = (R + R') / 2;
aRef = aRef(:);
switch method
    case "mvdr"
        [w, diagnostics] = solveClosedForm(R, aRef, cfg, "mvdr");
    case "dlmvdr"
        [w, diagnostics] = solveClosedForm(R, aRef, cfg, "dlmvdr");
    case "robust"
        [w, diagnostics] = solveRobust(R, aRef, cfg);
    otherwise
        error('uav:beamforming:UnsupportedMethod', ...
            'Supported methods are "mvdr", "dlmvdr", and "robust".');
end
end

function [w, diagnostics] = solveClosedForm(R, aRef, cfg, method)
M = size(R, 1);
if method == "mvdr"
    loading = minimumNumericalLoading(R, cfg);
else
    loading = cfg.beamforming.diagonalLoading * covarianceScale(R);
end
loadedCovariance = R + loading * eye(M);
[w, normalizer] = normalizedWeight(loadedCovariance, aRef);
diagnostics = struct( ...
    'method', method, ...
    'loading', loading, ...
    'loadedCovariance', loadedCovariance, ...
    'reciprocalCondition', rcond(loadedCovariance), ...
    'normalizer', normalizer);
end

function [w, diagnostics] = solveRobust(R, aRef, cfg)
epsilonS = cfg.beamforming.epsilonS;
gamma = cfg.beamforming.gamma;
referenceNorm = norm(aRef);
if epsilonS >= referenceNorm
    error('uav:beamforming:InfeasibleUncertainty', ...
        'Robust feasibility requires epsilonS < ||aRef||.');
end
if epsilonS == 0 && gamma == 0
    [w, diagnostics] = solveClosedForm(R, aRef, cfg, "mvdr");
    diagnostics.method = "robust";
    diagnostics.gamma = gamma;
    diagnostics.epsilonS = epsilonS;
    diagnostics.responseMargin = real(w' * aRef);
    diagnostics.usedMvdrRegressionPath = true;
    diagnostics.exitflag = 1;
    diagnostics.realSolution = [real(w); imag(w)];
    return;
end

M = size(R, 1);
regularizedCovariance = R + gamma * eye(M);
numericalLoading = minimumNumericalLoading(regularizedCovariance, cfg);
regularizedCovariance = regularizedCovariance + numericalLoading * eye(M);
Q = [real(regularizedCovariance), -imag(regularizedCovariance); ...
    imag(regularizedCovariance), real(regularizedCovariance)];
Q = (Q + Q') / 2;
realReference = [real(aRef); imag(aRef)];
solver = cfg.beamforming.solver;
interiorScale = (1 + solver.initialInteriorMargin) / ...
    (referenceNorm * (referenceNorm - epsilonS));
initial = interiorScale * realReference;
options = optimoptions('fmincon', ...
    'Algorithm', 'interior-point', ...
    'Display', 'off', ...
    'SpecifyObjectiveGradient', true, ...
    'SpecifyConstraintGradient', true, ...
    'OptimalityTolerance', solver.optimalityTolerance, ...
    'ConstraintTolerance', solver.constraintTolerance, ...
    'StepTolerance', solver.stepTolerance, ...
    'MaxIterations', solver.maxIterations);
[realSolution, objectiveValue, exitflag, output] = fmincon( ...
    @(x) realObjective(x, Q), initial, [], [], [], [], [], [], ...
    @(x) robustConstraint(x, realReference, epsilonS), options);
w = realSolution(1:M) + 1i * realSolution(M + 1:end);
responseMargin = real(w' * aRef) - epsilonS * norm(w);
if exitflag <= 0
    error('uav:beamforming:SolverFailure', ...
        'Robust optimization failed with exit flag %d.', exitflag);
end
if any(~isfinite(w)) || ~isfinite(responseMargin)
    error('uav:beamforming:NonfiniteWeight', ...
        'Robust optimization returned nonfinite weights.');
end
if responseMargin < 1 - solver.feasibilityTolerance
    error('uav:beamforming:ConstraintViolation', ...
        'Robust response margin %.16g is infeasible.', responseMargin);
end
diagnostics = struct( ...
    'method', "robust", ...
    'gamma', gamma, ...
    'epsilonS', epsilonS, ...
    'responseMargin', responseMargin, ...
    'usedMvdrRegressionPath', false, ...
    'exitflag', exitflag, ...
    'output', output, ...
    'objectiveValue', objectiveValue, ...
    'realSolution', realSolution, ...
    'numericalLoading', numericalLoading, ...
    'regularizedCovariance', regularizedCovariance);
end

function [objective, gradient] = realObjective(x, Q)
objective = x' * Q * x;
gradient = 2 * Q * x;
end

function [inequality, equality, inequalityGradient, equalityGradient] = ...
        robustConstraint(x, realReference, epsilonS)
vectorNorm = norm(x);
inequality = 1 - realReference' * x + epsilonS * vectorNorm;
equality = [];
if vectorNorm > 0
    inequalityGradient = -realReference + epsilonS * x / vectorNorm;
else
    inequalityGradient = -realReference;
end
equalityGradient = [];
end

function [w, normalizer] = normalizedWeight(covariance, aRef)
rawWeight = covariance \ aRef;
if any(~isfinite(rawWeight)) || norm(rawWeight) == 0
    rawWeight = pinv(covariance) * aRef;
end
normalizer = aRef' * rawWeight;
if ~isfinite(normalizer) || abs(normalizer) <= ...
        eps(max(1, norm(aRef) * norm(rawWeight)))
    error('uav:beamforming:DegenerateReference', ...
        'The reference steering vector cannot be normalized against R.');
end
w = rawWeight / normalizer;
end

function validateInputs(R, aRef, cfg)
if ~ismatrix(R) || size(R, 1) ~= size(R, 2) || any(~isfinite(R), 'all')
    error('uav:beamforming:InvalidCovariance', ...
        'R must be a finite square covariance matrix.');
end
if ~isvector(aRef) || numel(aRef) ~= size(R, 1) || any(~isfinite(aRef))
    error('uav:beamforming:InvalidReference', ...
        'aRef must be a finite vector matching R.');
end
if ~isstruct(cfg) || ~isfield(cfg, 'beamforming') || ...
        ~isfield(cfg.beamforming, 'minimumReciprocalCondition')
    error('uav:beamforming:InvalidConfiguration', ...
        'cfg is missing beamforming fields.');
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
