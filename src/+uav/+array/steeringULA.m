function a = steeringULA(thetaRad, sensorCount, spacingOverLambda)
%STEERINGULA 生成均匀线阵的窄带远场导向矢量。
%   thetaRad         - 1×Q或Q×1入射角，阵列法线方向为0；
%   sensorCount      - 阵元数M；
%   spacingOverLambda- 阵元间距与波长之比；
%   a                - M×Q导向矢量矩阵，每列范数为sqrt(M)。
%
%   相位约定为exp(j*2*pi*d/lambda*m*sin(theta))，与原论文半波长
%   ULA仿真及本工程w^H*x接收约定一致。

if ~isvector(thetaRad) || isempty(thetaRad) || any(~isfinite(thetaRad)) || ...
        any(abs(real(thetaRad)) > pi / 2 + 1e-12)
    error('uav:array:InvalidUlaAngles', ...
        'thetaRad must contain finite angles in [-pi/2, pi/2].');
end
if ~isscalar(sensorCount) || sensorCount < 1 || ...
        sensorCount ~= floor(sensorCount)
    error('uav:array:InvalidUlaSensorCount', ...
        'sensorCount must be a positive integer.');
end
if ~isscalar(spacingOverLambda) || ~isfinite(spacingOverLambda) || ...
        spacingOverLambda <= 0
    error('uav:array:InvalidUlaSpacing', ...
        'spacingOverLambda must be positive and finite.');
end

thetaRad = reshape(real(thetaRad), 1, []);
sensorIndices = (0:sensorCount - 1).';
a = exp(1i * 2 * pi * spacingOverLambda * ...
    sensorIndices .* sin(thetaRad));
end
