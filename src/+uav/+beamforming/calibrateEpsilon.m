function [epsilonS, diagnostics] = calibrateEpsilon(cfg, stream)
%CALIBRATEEPSILON 离线标定球形导向矢量失配半径epsilonS。
%   本函数只使用配置给出的“设计抖动分布”，不读取任何当前试次真值。
%   对每个设计姿态样本计算||a(eta)-a_ref||_2，再取配置分位数作为
%   epsilonS。这样同一实验配置中的所有方法、扫描点和试次共享同一半径。

if nargin < 2 || isempty(stream)
    stream = uav.util.makeStream(cfg.random.baseSeed, ...
        cfg.beamforming.calibration.substream);
end
validateInputs(cfg, stream);
calibration = cfg.beamforming.calibration;
sampleCount = calibration.sampleCount;
% 三轴设计标准差只用于鲁棒半径标定，不会覆盖实际实验抖动参数。
designStdRad = deg2rad(reshape(calibration.designJitterStdDeg, 3, 1));
offsets = designStdRad .* randn(stream, 3, sampleCount);
nominalAttitudeRad = deg2rad(reshape(cfg.attitude.nominalDeg, 3, 1));
desiredWorldDirection = cfg.derived.worldUnitVectors(:, 1);
aRef = cfg.derived.nominalSteeringReference;
mismatchNorms = zeros(sampleCount, 1);
for index = 1:sampleCount
    % 对期望信号执行与主场景相同的姿态->方向->流形映射。
    rotation = uav.geometry.rotationZYX(nominalAttitudeRad + offsets(:, index));
    arrayDirection = uav.geometry.worldToArray(rotation, desiredWorldDirection);
    directionCosines = uav.geometry.directionCoordinates(arrayDirection);
    steering = uav.array.steeringUPA(directionCosines, cfg.array.Nx, ...
        cfg.array.Ny, cfg.array.spacingOverLambda);
    mismatchNorms(index) = norm(steering - aRef);
end
% 使用经验分位数而不是当前试次最大误差，避免泄漏测试真值。
epsilonS = quantile(mismatchNorms, calibration.quantile);
if ~isfinite(epsilonS) || epsilonS < 0 || epsilonS >= norm(aRef)
    error('uav:beamforming:InfeasibleCalibration', ...
        'The calibrated epsilon must be finite and smaller than ||aRef||.');
end
diagnostics = struct( ...
    'designJitterStdDeg', calibration.designJitterStdDeg, ...
    'quantile', calibration.quantile, ...
    'sampleCount', sampleCount, ...
    'configuredSubstream', calibration.substream, ...
    'mismatchNorms', mismatchNorms, ...
    'referenceNorm', norm(aRef));
end

function validateInputs(cfg, stream)
if ~isstruct(cfg) || ~isfield(cfg, 'beamforming') || ...
        ~isfield(cfg.beamforming, 'calibration') || ~isfield(cfg, 'derived') || ...
        ~isfield(cfg, 'attitude') || ~isfield(cfg, 'array')
    error('uav:beamforming:InvalidConfiguration', ...
        'cfg is missing offline calibration fields.');
end
if ~isa(stream, 'RandStream')
    error('uav:beamforming:InvalidStream', ...
        'Calibration requires an explicit RandStream.');
end
end
