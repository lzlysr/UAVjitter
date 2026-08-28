function scenario = generate(cfg, stream)
%GENERATE 生成整段姿态轨迹和逐周期真实阵列流形。
%   cfg      - 已校验配置；
%   stream   - 显式RandStream，保证不同试次可以确定性复现；
%   scenario - 只保存姿态、旋转和方向等较小状态，不预生成M×L×K快拍。
%
%   姿态偏差采用三轴独立一阶高斯—马尔可夫过程：
%     delta_eta(k)=F*delta_eta(k-1)+nu(k)
%     Qnu=Qeta-F*Qeta*F'
%   初值直接从稳态分布N(0,Qeta)采样，因此无需额外烧入周期。

validateStream(stream);
K = cfg.timeline.K;
M = cfg.derived.M;
sourceCount = cfg.derived.J + 1;
F = diag(cfg.attitude.rho(:));
stdRad = deg2rad(cfg.attitude.stdDeg(:));
Qeta = diag(stdRad .^ 2);
% 用稳态协方差关系反推创新噪声协方差。
Qnu = Qeta - F * Qeta * F.';
if any(diag(Qnu) < -eps(max(1, max(diag(Qeta)))))
    error('uav:scenario:InvalidAttitudeModel', 'Stationary innovation covariance is not PSD.');
end
innovationStd = sqrt(max(diag(Qnu), 0));

attitudeOffsetRad = zeros(3, K);
% 第一个时刻从目标稳态分布采样，后续按AR(1)状态方程递推。
attitudeOffsetRad(:, 1) = stdRad .* randn(stream, 3, 1);
for timeIndex = 2:K
    attitudeOffsetRad(:, timeIndex) = F * attitudeOffsetRad(:, timeIndex - 1) + ...
        innovationStd .* randn(stream, 3, 1);
end

nominalRad = deg2rad(cfg.attitude.nominalDeg(:));
attitudeRad = nominalRad + attitudeOffsetRad;
% 以下数组按“坐标维×信号源×时间”组织，最后一维始终是周期索引。
rotations = zeros(3, 3, K);
arrayUnitVectors = zeros(3, sourceCount, K);
directionCosines = zeros(2, sourceCount, K);
azElRad = zeros(2, sourceCount, K);
steeringTruth = complex(zeros(M, sourceCount, K));
for timeIndex = 1:K
    % 严格执行 eta -> R -> u_array -> xi -> a 的物理映射链。
    R = uav.geometry.rotationZYX(attitudeRad(:, timeIndex));
    uArray = uav.geometry.worldToArray(R, cfg.derived.worldUnitVectors);
    [xi, directions] = uav.geometry.directionCoordinates(uArray);

    rotations(:, :, timeIndex) = R;
    arrayUnitVectors(:, :, timeIndex) = uArray;
    directionCosines(:, :, timeIndex) = xi;
    azElRad(:, :, timeIndex) = directions;
    steeringTruth(:, :, timeIndex) = uav.array.steeringUPA(xi, cfg.array.Nx, ...
        cfg.array.Ny, cfg.array.spacingOverLambda);
end

% 同时保留状态矩阵和协方差，便于模型验证实验检查。
scenario = struct( ...
    'F', F, ...
    'Qeta', Qeta, ...
    'Qnu', Qnu, ...
    'attitudeOffsetRad', attitudeOffsetRad, ...
    'attitudeRad', attitudeRad, ...
    'rotations', rotations, ...
    'arrayUnitVectors', arrayUnitVectors, ...
    'directionCosines', directionCosines, ...
    'azElRad', azElRad, ...
    'steeringTruth', steeringTruth);
end

function validateStream(stream)
if ~isa(stream, 'RandStream')
    error('uav:scenario:MissingStream', 'A RandStream input is required.');
end
end
