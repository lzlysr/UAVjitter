function [record, nextState] = runMethods(obs, priorState, cfg)
%RUNMETHODS Evaluate the fixed six methods on one shared observation.

methodNames = reshape(cfg.experiments.methodNames, 1, []);
if ~isequal(methodNames, ...
        ["oracle", "frozen", "scm", "dl", "static", "temporal"])
    error('uav:experiment:InvalidMethodContract', ...
        'The six-method runner requires the fixed configured method order.');
end
if ~isstruct(obs) || ~all(isfield(obs, {'X', 'Rx', 'truth'}))
    error('uav:experiment:InvalidObservation', ...
        'The fair runner requires one generated observation and isolated truth.');
end
if isempty(priorState)
    priorState = struct('frozenWeight', [], 'temporalState', []);
elseif ~isstruct(priorState) || ...
        ~all(isfield(priorState, {'frozenWeight', 'temporalState'}))
    error('uav:experiment:InvalidRunnerState', ...
        'Runner state must preserve frozen and temporal state separately.');
end

methodCount = numel(methodNames);
sinrLinear = zeros(1, methodCount);
desiredGain = zeros(1, methodCount);
interferenceLeakage = zeros(1, methodCount);
runtimeSec = zeros(1, methodCount);
incmNmse = nan(1, methodCount);
incmMinimumEigenvalue = nan(1, methodCount);
incmConditionNumber = nan(1, methodCount);
weights = complex(zeros(cfg.derived.M, methodCount));
diagnostics = cell(1, methodCount);
aRef = cfg.derived.nominalSteeringReference;
algorithmObservation = uav.experiment.algorithmObservation(obs);

for methodIndex = 1:methodCount
    methodTimer = tic;
    method = methodNames(methodIndex);
    estimatedIncm = [];
    switch method
        case "oracle"
            [weight, beamformingDiagnostics] = uav.beamforming.solve( ...
                "mvdr", obs.truth.trueIncm, obs.truth.trueASVs(:, 1), cfg);
            methodDiagnostics = struct('beamforming', beamformingDiagnostics, ...
                'truthAccess', "oracle-only");
        case "frozen"
            if isempty(priorState.frozenWeight)
                [priorState.frozenWeight, beamformingDiagnostics] = ...
                    uav.beamforming.solve("mvdr", obs.Rx, aRef, cfg);
                initialized = true;
            else
                beamformingDiagnostics = struct('method', "mvdr", ...
                    'source', "first-period-frozen");
                initialized = false;
            end
            weight = priorState.frozenWeight;
            methodDiagnostics = struct('beamforming', beamformingDiagnostics, ...
                'initializedThisPeriod', initialized);
        case "scm"
            [weight, beamformingDiagnostics] = uav.beamforming.solve( ...
                "mvdr", obs.Rx, aRef, cfg);
            methodDiagnostics = struct('beamforming', beamformingDiagnostics);
        case "dl"
            [weight, beamformingDiagnostics] = uav.beamforming.solve( ...
                "dlmvdr", obs.Rx, aRef, cfg);
            methodDiagnostics = struct('beamforming', beamformingDiagnostics);
        case "static"
            [estimatedIncm, ~, reconstructionDiagnostics] = ...
                uav.incm.reconstruct("static", algorithmObservation, [], cfg);
            [weight, beamformingDiagnostics] = uav.beamforming.solve( ...
                "robust", estimatedIncm, aRef, cfg);
            methodDiagnostics = struct('reconstruction', reconstructionDiagnostics, ...
                'beamforming', beamformingDiagnostics);
        case "temporal"
            [estimatedIncm, priorState.temporalState, reconstructionDiagnostics] = ...
                uav.incm.reconstruct("temporal", algorithmObservation, ...
                priorState.temporalState, cfg);
            [weight, beamformingDiagnostics] = uav.beamforming.solve( ...
                "robust", estimatedIncm, aRef, cfg);
            methodDiagnostics = struct('reconstruction', reconstructionDiagnostics, ...
                'beamforming', beamformingDiagnostics);
    end
    runtimeSec(methodIndex) = toc(methodTimer);
    metricTruth = obs.truth;
    if ~isempty(estimatedIncm)
        metricTruth.estimatedIncm = estimatedIncm;
    end
    metrics = uav.metrics.evaluate(weight, metricTruth, cfg);
    sinrLinear(methodIndex) = metrics.sinrLinear;
    desiredGain(methodIndex) = metrics.desiredGain;
    interferenceLeakage(methodIndex) = metrics.interferenceLeakage;
    if ~isempty(estimatedIncm)
        incmNmse(methodIndex) = metrics.incmNmse;
        eigenvalues = real(eig((estimatedIncm + estimatedIncm') / 2));
        incmMinimumEigenvalue(methodIndex) = min(eigenvalues);
        incmConditionNumber(methodIndex) = cond(estimatedIncm);
    end
    weights(:, methodIndex) = weight;
    diagnostics{methodIndex} = methodDiagnostics;
end

signature = uav.util.snapshotSignature(obs.X);
record = struct( ...
    'methodNames', methodNames, ...
    'sinrLinear', sinrLinear, ...
    'desiredGain', desiredGain, ...
    'interferenceLeakage', interferenceLeakage, ...
    'runtimeSec', runtimeSec, ...
    'incmNmse', incmNmse, ...
    'incmMinimumEigenvalue', incmMinimumEigenvalue, ...
    'incmConditionNumber', incmConditionNumber, ...
    'weights', weights, ...
    'diagnostics', {diagnostics}, ...
    'snapshotSignatures', repmat(signature, methodCount, 1), ...
    'epsilonS', cfg.beamforming.epsilonS, ...
    'gamma', cfg.beamforming.gamma);
nextState = priorState;
end
