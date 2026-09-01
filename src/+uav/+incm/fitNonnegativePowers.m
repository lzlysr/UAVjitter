function [powers, diagnostics] = fitNonnegativePowers( ...
        Rx, candidateManifold, noisePower, lambdaFit)
%FITNONNEGATIVEPOWERS 以非负Hermitian协方差拟合估计候选方向功率。
%   每个字典列对应a_q*a_q^H。矩阵经保Frobenius范数的实值
%   Hermitian向量化后，用带L2正则项的非负最小二乘估计候选功率。
%   本函数只使用样本协方差、候选流形和配置噪声，不读取任何真值。

validateInputs(Rx, candidateManifold, noisePower, lambdaFit);
M = size(Rx, 1);
Q = size(candidateManifold, 2);

inputHermitianResidual = norm(Rx - Rx', 'fro');
inputScale = max(norm(Rx, 'fro'), eps);
inputWasHermitian = inputHermitianResidual <= sqrt(eps) * inputScale;
residualCovariance = (Rx + Rx') / 2 - noisePower * eye(M);
target = hermitianVectorize(residualCovariance);

dictionary = zeros(numel(target), Q);
for candidateIndex = 1:Q
    candidate = candidateManifold(:, candidateIndex);
    dictionary(:, candidateIndex) = hermitianVectorize(candidate * candidate');
end

% 增广形式同时实现L2正则化和候选功率非负约束。
augmentedDictionary = [dictionary; sqrt(lambdaFit) * eye(Q)];
augmentedTarget = [target; zeros(Q, 1)];
powers = lsqnonneg(augmentedDictionary, augmentedTarget);

fitResidual = target - dictionary * powers;
regularizedResidual = augmentedTarget - augmentedDictionary * powers;
diagnostics = struct( ...
    'residualNorm', norm(fitResidual), ...
    'relativeResidual', norm(fitResidual) / max(norm(target), eps), ...
    'regularizedResidualNorm', norm(regularizedResidual), ...
    'dictionaryRank', rank(dictionary), ...
    'inputHermitianResidual', inputHermitianResidual, ...
    'inputWasHermitian', inputWasHermitian, ...
    'inputHermitianized', ~inputWasHermitian, ...
    'usedTruth', false);
end

function vector = hermitianVectorize(covariance)
%HERMITIANVECTORIZE 将Hermitian矩阵映射到保Frobenius范数的实向量。
M = size(covariance, 1);
lowerMask = tril(true(M), -1);
vector = [real(diag(covariance)); ...
    sqrt(2) * real(covariance(lowerMask)); ...
    sqrt(2) * imag(covariance(lowerMask))];
end

function validateInputs(Rx, candidateManifold, noisePower, lambdaFit)
%VALIDATEINPUTS 拒绝破坏协方差拟合数学合同的输入。
if ~isnumeric(Rx) || ~ismatrix(Rx) || size(Rx, 1) ~= size(Rx, 2) || ...
        isempty(Rx) || any(~isfinite(Rx), 'all')
    error('uav:incm:InvalidCovariance', ...
        'Rx must be a finite nonempty square covariance matrix.');
end
M = size(Rx, 1);
if ~isnumeric(candidateManifold) || ~ismatrix(candidateManifold) || ...
        size(candidateManifold, 1) ~= M || size(candidateManifold, 2) < 1 || ...
        any(~isfinite(candidateManifold), 'all')
    error('uav:incm:InvalidCandidateManifold', ...
        'candidateManifold must be a finite M-by-Q matrix with Q at least one.');
end
if ~isnumeric(noisePower) || ~isreal(noisePower) || ~isscalar(noisePower) || ...
        ~isfinite(noisePower) || noisePower < 0
    error('uav:incm:InvalidNoisePower', ...
        'noisePower must be a finite nonnegative real scalar.');
end
if ~isnumeric(lambdaFit) || ~isreal(lambdaFit) || ~isscalar(lambdaFit) || ...
        ~isfinite(lambdaFit) || lambdaFit < 0
    error('uav:incm:InvalidFitRegularization', ...
        'lambdaFit must be a finite nonnegative real scalar.');
end
end
