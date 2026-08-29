function [record, nextState] = runMethods(obs, priorState, cfg, methodNames)
%RUNMETHODS 在同一个观测上公平运行指定的主线比较方法。
%   priorState只包含冻结权值和时序INCM上一周期状态，二者相互独立。
%   record返回每种方法的指标、运行时间、权值和诊断；所有方法使用相同
%   aRef、gamma、epsilonS以及同一快拍X。

if nargin < 4
    % 保留原三参数接口，实验3、4仍按固定六方法合同运行。
    methodNames = cfg.experiments.methodNames;
end
methodNames = reshape(string(methodNames), 1, []);
allowedMethods = ["oracle", "frozen", "scm", "dl", "static", "temporal"];
if isempty(methodNames) || any(strlength(methodNames) == 0) || ...
        numel(unique(methodNames)) ~= numel(methodNames) || ...
        any(~ismember(methodNames, allowedMethods))
    error('uav:experiment:InvalidMethodSelection', ...
        'The requested methods must form a nonempty unique supported subset.');
end
if ~isstruct(obs) || ~all(isfield(obs, {'X', 'Rx', 'truth'}))
    error('uav:experiment:InvalidObservation', ...
        'The fair runner requires one generated observation and isolated truth.');
end
if isempty(priorState)
    % 首周期同时初始化冻结基线和时序方法；之后沿各自状态递推。
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
            % 唯一允许在求权值阶段读取真实INCM和真实期望流形的方法。
            [weight, beamformingDiagnostics] = uav.beamforming.solve( ...
                "mvdr", obs.truth.trueIncm, obs.truth.trueASVs(:, 1), cfg);
            methodDiagnostics = struct('beamforming', beamformingDiagnostics, ...
                'truthAccess', "oracle-only");
        case "frozen"
            % 首周期用SCM求一次权值，后续不随姿态和观测更新。
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
            % 直接把含期望信号的当前SCM代入标准MVDR。
            [weight, beamformingDiagnostics] = uav.beamforming.solve( ...
                "mvdr", obs.Rx, aRef, cfg);
            methodDiagnostics = struct('beamforming', beamformingDiagnostics);
        case "dl"
            % 当前SCM加显式尺度化对角加载。
            [weight, beamformingDiagnostics] = uav.beamforming.solve( ...
                "dlmvdr", obs.Rx, aRef, cfg);
            methodDiagnostics = struct('beamforming', beamformingDiagnostics);
        case "static"
            % 当前周期静态Capon-INCM + 与时序方法相同的鲁棒权值求解器。
            [estimatedIncm, ~, reconstructionDiagnostics] = ...
                uav.incm.reconstruct("static", algorithmObservation, [], cfg);
            [weight, beamformingDiagnostics] = uav.beamforming.solve( ...
                "robust", estimatedIncm, aRef, cfg);
            methodDiagnostics = struct('reconstruction', reconstructionDiagnostics, ...
                'beamforming', beamformingDiagnostics);
        case "temporal"
            % 使用并更新上一周期增强谱状态。
            [estimatedIncm, priorState.temporalState, reconstructionDiagnostics] = ...
                uav.incm.reconstruct("temporal", algorithmObservation, ...
                priorState.temporalState, cfg);
            [weight, beamformingDiagnostics] = uav.beamforming.solve( ...
                "robust", estimatedIncm, aRef, cfg);
            methodDiagnostics = struct('reconstruction', reconstructionDiagnostics, ...
                'beamforming', beamformingDiagnostics);
    end
    runtimeSec(methodIndex) = toc(methodTimer);
    % 真值只在权值已经确定后交给统一指标函数。
    metricTruth = obs.truth;
    if ~isempty(estimatedIncm)
        % 只有真正重构INCM的方法才报告NMSE、最小特征值和条件数。
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
% 同一摘要复制到全部所选方法，用于实验层自动检查配对观测。
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
