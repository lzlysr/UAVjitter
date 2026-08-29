function [w, diagnostics] = solve(method, R, aRef, cfg)
%SOLVE 统一求解MVDR、对角加载MVDR和鲁棒波束形成权值。
%   R    - M×M协方差或重构INCM；
%   aRef - M×1名义期望导向矢量；
%   w    - M×1复权值，采用w^H*x的接收约定。
%
%   robust方法求解：
%     min_w w^H*R*w + gamma*||w||_2^2
%     s.t. Re{w^H*aRef}-epsilonS*||w||_2 >= 1。

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
    % 标准MVDR只在数值病态时加入满足条件数要求的最小加载。
    loading = minimumNumericalLoading(R, cfg);
else
    % DL-MVDR使用与协方差平均功率成比例的显式加载。
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
inputScale = max(norm(R, 2), 1);
psdTolerance = cfg.beamforming.covariancePsdTolerance * inputScale;
minimumInputEigenvalue = min(real(eig(R)));
if minimumInputEigenvalue < -psdTolerance
    error('uav:beamforming:IndefiniteCovariance', ...
        'Robust R must be PSD within the configured relative tolerance.');
end
if epsilonS >= referenceNorm
    error('uav:beamforming:InfeasibleUncertainty', ...
        'Robust feasibility requires epsilonS < ||aRef||.');
end
if epsilonS == 0
    % 合同要求：不确定集半径为零时必须严格返回标准MVDR，gamma不参与。
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
% gamma*||w||^2等价于在目标协方差上增加gamma*I。
objectiveCovariance = R + gamma * eye(M);
% 将复Hermitian二次型转换成2M维实对称二次型。
Q = [real(objectiveCovariance), -imag(objectiveCovariance); ...
    imag(objectiveCovariance), real(objectiveCovariance)];
Q = (Q + Q') / 2;
realReference = [real(aRef); imag(aRef)];
solver = cfg.beamforming.solver;
% 沿aRef方向构造严格位于可行域内部的初始点，帮助内点法收敛。
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
% 用原始复数形式重新计算约束余量，不能只依赖优化器退出标志。
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
    'minimumInputEigenvalue', minimumInputEigenvalue, ...
    'psdTolerance', psdTolerance, ...
    'objectiveCovariance', objectiveCovariance, ...
    'realQuadraticMatrix', Q);
end

function [objective, gradient] = realObjective(x, Q)
% 实变量二次目标及其解析梯度，供fmincon避免数值差分。
objective = x' * Q * x;
gradient = 2 * Q * x;
end

function [inequality, equality, inequalityGradient, equalityGradient] = ...
        robustConstraint(x, realReference, epsilonS)
% fmincon使用c(x)<=0，因此原约束改写为1-a^T*x+epsilon*||x||<=0。
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
% 闭式MVDR：w=R^{-1}a/(a^H R^{-1}a)，保证w^H*a=1。
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
% 计算协方差矩阵的平均功率，避免在病态情况下使用最大特征值。
scale = max(real(trace(R)) / size(R, 1), eps);
end

function loading = minimumNumericalLoading(R, cfg)
% 在不必要时保持loading=0；病态时结合特征值界和1范数界估计加载量。
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
