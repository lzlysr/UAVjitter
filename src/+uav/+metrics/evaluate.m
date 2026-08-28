function metrics = evaluate(w, truth, cfg)
%EVALUATE 使用隔离真值统一计算波束形成性能指标。
%   本函数是非Oracle算法可以间接使用真值的唯一评价位置之一：算法权值
%   已经求出后，才把w与真实导向矢量、真实INCM结合。所有功率先在线性域
%   计算；dB值只作为输出字段，蒙特卡洛统计仍使用sinrLinear。

validateInputs(w, truth, cfg);
desiredAsv = truth.trueASVs(:, 1);
desiredGain = abs(w' * desiredAsv)^2;
desiredOutputPower = truth.desiredPower * desiredGain;
% 分母包含全部干扰和接收机噪声。
interferencePlusNoisePower = real(w' * truth.trueIncm * w);
if interferencePlusNoisePower <= 0
    error('uav:metrics:InvalidIncm', ...
        'trueIncm must yield positive interference-plus-noise output power.');
end

if isempty(truth.interfererPowers)
    % J=0时“干扰泄漏”定义为0，但SINR分母仍包含噪声。
    interferenceLeakage = 0;
else
    interfererAsvs = truth.trueASVs(:, 2:end);
    interferenceLeakage = real(w' * interfererAsvs * ...
        diag(truth.interfererPowers(:)) * interfererAsvs' * w);
end

metrics = struct( ...
    'sinrLinear', desiredOutputPower / interferencePlusNoisePower, ...
    'sinrDb', 10 * log10(desiredOutputPower / interferencePlusNoisePower), ...
    'desiredGain', desiredGain, ...
    'interferenceLeakage', interferenceLeakage, ...
    'desiredOutputPower', desiredOutputPower, ...
    'interferencePlusNoisePower', interferencePlusNoisePower);
if isfield(truth, 'estimatedIncm')
    % 只有INCM方法传入估计矩阵时才附加归一化Frobenius误差。
    estimate = truth.estimatedIncm;
    if ~isequal(size(estimate), size(truth.trueIncm)) || any(~isfinite(estimate), 'all')
        error('uav:metrics:InvalidIncmEstimate', ...
            'estimatedIncm must be finite and match trueIncm.');
    end
    metrics.incmNmse = norm(estimate - truth.trueIncm, 'fro')^2 / ...
        norm(truth.trueIncm, 'fro')^2;
end
end

function validateInputs(w, truth, cfg)
if ~isstruct(cfg)
    error('uav:metrics:InvalidConfiguration', 'cfg must be a configuration struct.');
end
required = {'trueASVs', 'desiredPower', 'interfererPowers', 'trueIncm'};
if ~isstruct(truth) || ~all(isfield(truth, required))
    error('uav:metrics:InvalidTruth', 'truth does not contain the required metric fields.');
end
M = size(truth.trueIncm, 1);
if ~isequal(size(truth.trueIncm), [M, M]) || size(truth.trueASVs, 1) ~= M || ...
        ~isvector(w) || numel(w) ~= M || any(~isfinite(w)) || ...
        any(~isfinite(truth.trueASVs), 'all') || any(~isfinite(truth.trueIncm), 'all') || ...
        ~isscalar(truth.desiredPower) || ~isfinite(truth.desiredPower) || truth.desiredPower < 0 || ...
        numel(truth.interfererPowers) ~= size(truth.trueASVs, 2) - 1 || ...
        any(~isfinite(truth.interfererPowers)) || any(truth.interfererPowers < 0)
    error('uav:metrics:InvalidTruth', 'Metric inputs have incompatible dimensions or values.');
end
end
