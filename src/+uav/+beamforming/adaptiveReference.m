function [w, nextState, diagnostics] = adaptiveReference(method, obs, priorState, cfg)
%ADAPTIVEREFERENCE 使用已知期望训练序列实现经典NLMS或RLS。
%   每个快拍的参考信号为sqrt(Ps)*pilot，使理想无失真波束输出与接收
%   期望信号幅度一致。priorState只保存上一周期权值及RLS逆相关矩阵，
%   不包含真实方向、真实协方差或干扰符号。

method = lower(string(method));
validateInputs(obs, cfg);
X = obs.X;
pilot = reshape(obs.pilot, 1, []);
desiredReference = sqrt(cfg.derived.desiredPower) * pilot;
M = size(X, 1);

switch method
    case "nlms"
        w = initialWeight(priorState, M, "nlms");
        stepSize = cfg.beamforming.classic.nlms.stepSize;
        regularization = cfg.beamforming.classic.nlms.regularization;
        for snapshotIndex = 1:size(X, 2)
            x = X(:, snapshotIndex);
            errorSample = desiredReference(snapshotIndex) - w' * x;
            w = w + stepSize * x * conj(errorSample) / ...
                (regularization + real(x' * x));
        end
        nextState = struct('method', "nlms", 'weight', w);
        diagnostics = struct('method', "nlms", 'stepSize', stepSize, ...
            'regularization', regularization, 'snapshotCount', size(X, 2));
    case "rls"
        [w, inverseCorrelation] = initialRlsState(priorState, M, cfg);
        forgettingFactor = cfg.beamforming.classic.rls.forgettingFactor;
        for snapshotIndex = 1:size(X, 2)
            x = X(:, snapshotIndex);
            denominator = forgettingFactor + real(x' * inverseCorrelation * x);
            gain = inverseCorrelation * x / denominator;
            errorSample = desiredReference(snapshotIndex) - w' * x;
            w = w + gain * conj(errorSample);
            inverseCorrelation = (inverseCorrelation - ...
                gain * x' * inverseCorrelation) / forgettingFactor;
            inverseCorrelation = (inverseCorrelation + inverseCorrelation') / 2;
        end
        nextState = struct('method', "rls", 'weight', w, ...
            'inverseCorrelation', inverseCorrelation);
        diagnostics = struct('method', "rls", ...
            'forgettingFactor', forgettingFactor, ...
            'snapshotCount', size(X, 2));
    otherwise
        error('uav:beamforming:UnsupportedAdaptiveMethod', ...
            'The adaptive reference method must be "nlms" or "rls".');
end

if any(~isfinite(w))
    error('uav:beamforming:NonfiniteWeight', ...
        'The classical adaptive method returned nonfinite weights.');
end
diagnostics.finalTrainingMse = mean(abs(desiredReference - w' * X).^2);
end

function w = initialWeight(priorState, M, method)
if isempty(priorState)
    w = complex(zeros(M, 1));
    return;
end
if ~isstruct(priorState) || ~isfield(priorState, 'method') || ...
        priorState.method ~= method || ~isfield(priorState, 'weight') || ...
        ~isequal(size(priorState.weight), [M, 1])
    error('uav:beamforming:InvalidAdaptiveState', ...
        'The adaptive prior state does not match the requested method.');
end
w = priorState.weight;
end

function [w, inverseCorrelation] = initialRlsState(priorState, M, cfg)
if isempty(priorState)
    w = complex(zeros(M, 1));
    inverseCorrelation = cfg.beamforming.classic.rls.initialInverseScale * eye(M);
    return;
end
w = initialWeight(priorState, M, "rls");
if ~isfield(priorState, 'inverseCorrelation') || ...
        ~isequal(size(priorState.inverseCorrelation), [M, M]) || ...
        any(~isfinite(priorState.inverseCorrelation), 'all')
    error('uav:beamforming:InvalidAdaptiveState', ...
        'The RLS prior state must contain a finite inverse-correlation matrix.');
end
inverseCorrelation = priorState.inverseCorrelation;
end

function validateInputs(obs, cfg)
if ~isstruct(obs) || ~all(isfield(obs, {'X', 'pilot'})) || ...
        ~ismatrix(obs.X) || size(obs.X, 1) ~= cfg.derived.M || ...
        size(obs.X, 2) ~= numel(obs.pilot) || isempty(obs.pilot) || ...
        any(~isfinite(obs.X), 'all') || any(~isfinite(obs.pilot))
    error('uav:beamforming:InvalidTrainingObservation', ...
        'NLMS/RLS require finite M-by-L snapshots and a matching pilot.');
end
end
