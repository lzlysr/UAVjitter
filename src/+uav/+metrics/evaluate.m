function metrics = evaluate(w, truth, cfg)
%EVALUATE Report output SINR and gain/leakage metrics against isolated truth.

validateInputs(w, truth, cfg);
desiredAsv = truth.trueASVs(:, 1);
desiredGain = abs(w' * desiredAsv)^2;
desiredOutputPower = truth.desiredPower * desiredGain;
interferencePlusNoisePower = real(w' * truth.trueIncm * w);
if interferencePlusNoisePower <= 0
    error('uav:metrics:InvalidIncm', ...
        'trueIncm must yield positive interference-plus-noise output power.');
end

if isempty(truth.interfererPowers)
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
