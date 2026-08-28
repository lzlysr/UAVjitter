function obs = generatePeriod(cfg, scenario, k, stream)
%GENERATEPERIOD 生成指定处理周期的窄带阵列快拍。
%   cfg/scenario - 配置与generate生成的动态场景；
%   k            - 当前周期索引；
%   stream       - 当前试次、当前周期专用随机流；
%   obs.X        - M×L复快拍矩阵；
%   obs.Rx       - M×M样本协方差X*X'/L；
%   obs.pilot    - 1×L已知期望QPSK训练序列；
%   obs.truth    - 真实阵列流形和真实INCM，只允许Oracle和指标模块读取。

validateStream(stream);
if ~isscalar(k) || k < 1 || k > cfg.timeline.K || k ~= floor(k)
    error('uav:scenario:InvalidPeriodIndex', 'k must index one configured period.');
end
if ~isfield(scenario, 'steeringTruth') || size(scenario.steeringTruth, 3) < k
    error('uav:scenario:InvalidScenario', 'scenario does not contain steering truth for k.');
end

M = cfg.derived.M;
J = cfg.derived.J;
L = cfg.timeline.L;
trueASVs = scenario.steeringTruth(:, :, k);
% 期望信号和每个邻区干扰使用相互独立的单位功率QPSK序列。
pilot = qpsk(stream, 1, L);
interfererSymbols = qpsk(stream, J, L);
% CN(0,sigma_n^2 I)可由独立实部/虚部N(0,sigma_n^2/2)构造。
noise = sqrt(cfg.derived.noisePower / 2) * (randn(stream, M, L) + ...
    1i * randn(stream, M, L));
desiredTerm = sqrt(cfg.derived.desiredPower) * trueASVs(:, 1) * pilot;
if J == 0
    % 空干扰场景仍保留噪声INCM，避免出现空矩阵或特殊指标定义。
    interferenceTerm = zeros(M, L);
    trueIncm = cfg.derived.noisePower * eye(M);
else
    interfererASVs = trueASVs(:, 2:end);
    weightedSymbols = sqrt(cfg.derived.interfererPowers(:)) .* interfererSymbols;
    interferenceTerm = interfererASVs * weightedSymbols;
    trueIncm = cfg.derived.noisePower * eye(M) + interfererASVs * ...
        diag(cfg.derived.interfererPowers) * interfererASVs';
end
% 窄带等效基带观测：期望项+同频干扰项+空间白噪声。
X = desiredTerm + interferenceTerm + noise;

% 真值集中封装，后续由experiment.algorithmObservation显式剥离。
truth = struct( ...
    'trueASVs', trueASVs, ...
    'trueIncm', trueIncm, ...
    'sourcePowers', cfg.derived.sourcePowers, ...
    'desiredPower', cfg.derived.desiredPower, ...
    'interfererPowers', cfg.derived.interfererPowers, ...
    'worldUnitVectors', cfg.derived.worldUnitVectors, ...
    'arrayUnitVectors', scenario.arrayUnitVectors(:, :, k), ...
    'directionCosines', scenario.directionCosines(:, :, k), ...
    'azElRad', scenario.azElRad(:, :, k), ...
    'interfererSymbols', interfererSymbols);
obs = struct( ...
    'X', X, ...
    'Rx', X * X' / L, ...
    'pilot', pilot, ...
    'k', k, ...
    'nominalSteeringReference', cfg.derived.nominalSteeringReference, ...
    'truth', truth);
end

function symbols = qpsk(stream, nRows, nColumns)
%QPSK 生成幅度为1、四个星座点等概率出现的独立QPSK符号。
symbols = ((2 * (rand(stream, nRows, nColumns) >= 0.5) - 1) + ...
    1i * (2 * (rand(stream, nRows, nColumns) >= 0.5) - 1)) / sqrt(2);
end

function validateStream(stream)
if ~isa(stream, 'RandStream')
    error('uav:scenario:MissingStream', 'A RandStream input is required.');
end
end
