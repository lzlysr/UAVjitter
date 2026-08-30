function [RinHat, diagnostics] = guLeshem2012Rec( ...
        sampleCovariance, scanManifold, quadratureWeights)
%GULESHEM2012REC 按Gu–Leshem 2012式(6)–(7)重构INCM。
%   sampleCovariance - M×M样本协方差\hat{R}；
%   scanManifold     - M×S角域补集扫描流形，每列为d(theta_s)；
%   quadratureWeights- 1×S或S×1数值求积权重；
%   RinHat           - sum_s weight_s*P_s*d_s*d_s^H。
%
%   本函数是原论文Rec方法的公式级核心，不实施对角加载、噪声底座扣除、
%   显式噪声单位阵、特征值修复或收缩。调用者必须提供可逆样本协方差；
%   若有限快拍导致矩阵奇异，应由独立迁移实验明确报告，不能静默改写原式。

validateInputs(sampleCovariance, scanManifold, quadratureWeights);
sampleCovariance = (sampleCovariance + sampleCovariance') / 2;
quadratureWeights = reshape(real(quadratureWeights), 1, []);

if rcond(sampleCovariance) <= eps(class(sampleCovariance))
    error('uav:incm:SingularSampleCovariance', ...
        'Gu-Leshem Rec requires an invertible sample covariance matrix.');
end

% 论文式(6)：逐列计算d^H*R^{-1}*d，不显式形成逆矩阵。
solvedManifold = sampleCovariance \ scanManifold;
denominator = real(sum(conj(scanManifold) .* solvedManifold, 1));
tolerance = eps(max(1, vecnorm(scanManifold, 2, 1) .* ...
    vecnorm(solvedManifold, 2, 1)));
if any(~isfinite(denominator)) || any(denominator <= tolerance)
    error('uav:incm:InvalidCaponDenominator', ...
        'Gu-Leshem Rec produced a nonpositive Capon denominator.');
end
caponPower = 1 ./ denominator;

% 论文式(7)的离散求积。权重只表示积分测度，不附加任何经验功率修正。
weightedPower = quadratureWeights .* caponPower;
RinHat = (scanManifold .* weightedPower) * scanManifold';
RinHat = (RinHat + RinHat') / 2;

diagnostics = struct( ...
    'caponPower', caponPower, ...
    'quadratureWeights', quadratureWeights, ...
    'weightedPower', weightedPower, ...
    'loading', 0, ...
    'usedPseudoInverse', false, ...
    'sampleCovarianceReciprocalCondition', rcond(sampleCovariance), ...
    'scanPointCount', size(scanManifold, 2));
end

function validateInputs(sampleCovariance, scanManifold, quadratureWeights)
%VALIDATEINPUTS 拒绝会改变原论文公式含义的缺失或非法输入。
if ~ismatrix(sampleCovariance) || ...
        size(sampleCovariance, 1) ~= size(sampleCovariance, 2) || ...
        isempty(sampleCovariance) || any(~isfinite(sampleCovariance), 'all')
    error('uav:incm:InvalidGuLeshemCovariance', ...
        'sampleCovariance must be a finite nonempty square matrix.');
end
M = size(sampleCovariance, 1);
if ~ismatrix(scanManifold) || size(scanManifold, 1) ~= M || ...
        size(scanManifold, 2) < 1 || any(~isfinite(scanManifold), 'all')
    error('uav:incm:InvalidGuLeshemManifold', ...
        'scanManifold must be a finite M-by-S matrix with S at least one.');
end
if ~isvector(quadratureWeights) || ...
        numel(quadratureWeights) ~= size(scanManifold, 2) || ...
        any(~isfinite(quadratureWeights)) || any(quadratureWeights < 0) || ...
        ~any(quadratureWeights > 0)
    error('uav:incm:InvalidGuLeshemWeights', ...
        'quadratureWeights must be finite, nonnegative, and match scan points.');
end
end
