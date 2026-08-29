function [w, diagnostics] = powerInversion(R, cfg)
%POWERINVERSION 实现固定参考阵元约束下的经典功率倒置。
%   求解min w^H R w，约束e_r^H w=1。该约束只固定一个参考阵元，
%   不保证期望信号方向无失真，因此能够体现PI可能压低期望信号的局限。

M = cfg.derived.M;
if ~isequal(size(R), [M, M]) || any(~isfinite(R), 'all')
    error('uav:beamforming:InvalidCovariance', ...
        'PI requires a finite M-by-M covariance matrix.');
end
R = (R + R') / 2;
referenceElement = cfg.beamforming.classic.pi.referenceElement;
referenceVector = zeros(M, 1);
referenceVector(referenceElement) = 1;
rawWeight = pinv(R) * referenceVector;
normalizer = referenceVector' * rawWeight;
if ~isfinite(normalizer) || abs(normalizer) <= eps(max(1, norm(rawWeight)))
    error('uav:beamforming:DegenerateReference', ...
        'The configured PI reference element cannot be normalized.');
end
w = rawWeight / normalizer;
diagnostics = struct('method', "pi", ...
    'referenceElement', referenceElement, ...
    'referenceResponse', referenceVector' * w, ...
    'covarianceRank', rank(R));
end
