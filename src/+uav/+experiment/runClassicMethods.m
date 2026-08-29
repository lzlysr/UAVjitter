function [record, nextState] = runClassicMethods(obs, priorState, cfg, methodNames)
%RUNCLASSICMETHODS 在同一观测上运行用户选择的经典方法组合。
%   methodNames可从配置支持的11种方法中任意选择且保持给定顺序。除
%   Oracle外，所有方法只读取白名单观测、配置先验和自身上一周期状态；
%   真实INCM和真实导向矢量只在Oracle求权值及统一指标评价阶段使用。

methodNames = reshape(string(methodNames), 1, []);
validateMethodSelection(methodNames);
if ~isstruct(obs) || ~all(isfield(obs, {'X', 'Rx', 'pilot', 'truth'}))
    error('uav:experiment:InvalidObservation', ...
        'Classic methods require one generated observation with isolated truth.');
end
state = initializeState(priorState);
safeObservation = uav.experiment.algorithmObservation(obs);
aRef = cfg.derived.nominalSteeringReference;
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
        case "nsb"
            [weight, beamformingDiagnostics] = ...
                uav.beamforming.nullSteering(cfg.derived.nominalSteering, cfg);
            methodDiagnostics = struct('beamforming', beamformingDiagnostics, ...
                'directionSource', "configured-nominal-directions");
        case "frozen_mvdr"
            if isempty(state.frozenWeight)
                [state.frozenWeight, beamformingDiagnostics] = ...
                    uav.beamforming.solve("mvdr", safeObservation.Rx, aRef, cfg);
                initialized = true;
            else
                beamformingDiagnostics = struct('method', "mvdr", ...
                    'source', "first-period-frozen");
                initialized = false;
            end
            weight = state.frozenWeight;
            methodDiagnostics = struct('beamforming', beamformingDiagnostics, ...
                'initializedThisPeriod', initialized);
        case "scm_mvdr"
            [weight, beamformingDiagnostics] = uav.beamforming.solve( ...
                "mvdr", safeObservation.Rx, aRef, cfg);
            methodDiagnostics = struct('beamforming', beamformingDiagnostics);
        case "nlms"
            [weight, state.nlms, beamformingDiagnostics] = ...
                uav.beamforming.adaptiveReference( ...
                "nlms", safeObservation, state.nlms, cfg);
            methodDiagnostics = struct('beamforming', beamformingDiagnostics, ...
                'referenceSource', "known-desired-pilot");
        case "rls"
            [weight, state.rls, beamformingDiagnostics] = ...
                uav.beamforming.adaptiveReference( ...
                "rls", safeObservation, state.rls, cfg);
            methodDiagnostics = struct('beamforming', beamformingDiagnostics, ...
                'referenceSource', "known-desired-pilot");
        case "pi"
            [weight, beamformingDiagnostics] = ...
                uav.beamforming.powerInversion(safeObservation.Rx, cfg);
            methodDiagnostics = struct('beamforming', beamformingDiagnostics);
        case "dl_mvdr"
            [weight, beamformingDiagnostics] = uav.beamforming.solve( ...
                "dlmvdr", safeObservation.Rx, aRef, cfg);
            methodDiagnostics = struct('beamforming', beamformingDiagnostics);
        case "esp"
            [weight, beamformingDiagnostics] = ...
                uav.beamforming.eigenspaceProjection( ...
                safeObservation.Rx, aRef, cfg);
            methodDiagnostics = struct('beamforming', beamformingDiagnostics);
        case "wco"
            [weight, beamformingDiagnostics] = uav.beamforming.solve( ...
                "wco", safeObservation.Rx, aRef, cfg);
            methodDiagnostics = struct('beamforming', beamformingDiagnostics);
        case "static_incm_mvdr"
            [estimatedIncm, ~, reconstructionDiagnostics] = ...
                uav.incm.reconstruct("static", safeObservation, [], cfg);
            reconstructionDiagnostics.method = "static";
            [weight, beamformingDiagnostics] = uav.beamforming.solve( ...
                "mvdr", estimatedIncm, aRef, cfg);
            methodDiagnostics = struct( ...
                'reconstruction', reconstructionDiagnostics, ...
                'beamforming', beamformingDiagnostics);
    end
    runtimeSec(methodIndex) = toc(methodTimer);
    metricTruth = obs.truth;
    if ~isempty(estimatedIncm)
        metricTruth.estimatedIncm = estimatedIncm;
    end
    methodMetrics = uav.metrics.evaluate(weight, metricTruth, cfg);
    sinrLinear(methodIndex) = methodMetrics.sinrLinear;
    desiredGain(methodIndex) = methodMetrics.desiredGain;
    interferenceLeakage(methodIndex) = methodMetrics.interferenceLeakage;
    if ~isempty(estimatedIncm)
        incmNmse(methodIndex) = methodMetrics.incmNmse;
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
    'snapshotSignatures', repmat(signature, methodCount, 1));
nextState = state;
end

function state = initializeState(priorState)
if isempty(priorState)
    state = struct('frozenWeight', [], 'nlms', [], 'rls', []);
    return;
end
required = {'frozenWeight', 'nlms', 'rls'};
if ~isstruct(priorState) || ~all(isfield(priorState, required))
    error('uav:experiment:InvalidClassicState', ...
        'Classic state must preserve frozen-MVDR, NLMS, and RLS states.');
end
state = priorState;
end

function validateMethodSelection(methodNames)
allowed = ["oracle", "nsb", "frozen_mvdr", "scm_mvdr", "nlms", ...
    "rls", "pi", "dl_mvdr", "esp", "wco", "static_incm_mvdr"];
if isempty(methodNames) || any(strlength(methodNames) == 0) || ...
        numel(unique(methodNames)) ~= numel(methodNames) || ...
        any(~ismember(methodNames, allowed))
    error('uav:experiment:InvalidClassicMethodSelection', ...
        'Classic methods must be a nonempty unique subset of supported names.');
end
end
