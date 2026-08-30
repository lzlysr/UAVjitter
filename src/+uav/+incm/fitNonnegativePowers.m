function [powers, diagnostics] = fitNonnegativePowers( ...
        Rx, candidateManifold, noisePower, lambdaFit)
%FITNONNEGATIVEPOWERS 以非负Hermitian协方差拟合估计候选方向功率。
%   Rx                - M×M有限样本协方差；允许含微小非Hermitian数值误差；
%   candidateManifold - M×Q有限候选流形，每列对应一个候选方向；
%   noisePower        - 非负有限标量噪声功率；
%   lambdaFit         - 非负有限标量L2正则化系数；
%   powers            - Q×1非负候选功率估计；
%   diagnostics       - 残差、字典秩和输入Hermitian化诊断。
%
%   每个字典列均为a_q*a_q^H的实值Hermitian向量化。实对角元素保留，
%   严格下三角的实部和虚部分别乘sqrt(2)，从而严格保持Frobenius范数。
%   本函数只使用传入协方差和候选流形，不读取任何真实量。

validateInputs(Rx, candidateManifold, noisePower, lambdaFit);
M = size(Rx, 1);
Q = size(candidateManifold, 2);

% 先数值Hermitian化，避免有限快拍计算中的舍入误差破坏实值拟合空间。
inputHermitianResidual = norm(Rx - Rx', 'fro');
inputScale = max(norm(Rx, 'fro'), eps);
inputWasHermitian = inputHermitianResidual <= sqrt(eps) * inputScale;
residualCovariance = (Rx + Rx') / 2 - noisePower * eye(M);
target = hermitianVectorize(residualCovariance);

dictionary = zeros(numel(target), Q);
for candidateIndex = 1:Q
    candidateCovariance = candidateManifold(:, candidateIndex) * ...
        candidateManifold(:, candidateIndex)';
    dictionary(:, candidateIndex) = hermitianVectorize(candidateCovariance);
end

% 将L2正则化写成增广NNLS，保持候选功率的非负物理约束。
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
%VALIDATEINPUTS 拒绝改变协方差拟合数学合同的无效输入。
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
