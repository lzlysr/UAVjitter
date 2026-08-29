function [w, diagnostics] = nullSteering(nominalSteering, cfg)
%NULLSTEERING 按标称方向实现经典零陷导向波束形成。
%   nominalSteering的第一列是期望信号导向矢量，其余列是干扰导向
%   矢量。权值满足期望方向单位响应、各标称干扰方向零响应，并在所有
%   可行权值中取二范数最小解。该方法不使用当前周期真实方向或INCM。

M = cfg.derived.M;
if ~ismatrix(nominalSteering) || size(nominalSteering, 1) ~= M || ...
        size(nominalSteering, 2) < 1 || any(~isfinite(nominalSteering), 'all')
    error('uav:beamforming:InvalidNominalSteering', ...
        'nominalSteering must be a finite M-by-sourceCount matrix.');
end

constraintCount = size(nominalSteering, 2);
desiredResponses = [1; zeros(constraintCount - 1, 1)];
gramMatrix = nominalSteering' * nominalSteering;
% pinv保留经典最小范数约束解，并允许约束矩阵接近奇异时稳定计算。
w = nominalSteering * pinv(gramMatrix) * desiredResponses;
actualResponses = nominalSteering' * w;
constraintResidual = norm(actualResponses - desiredResponses);
diagnostics = struct( ...
    'method', "nsb", ...
    'desiredResponses', desiredResponses, ...
    'actualResponses', actualResponses, ...
    'constraintResidual', constraintResidual, ...
    'constraintRank', rank(nominalSteering));
end
