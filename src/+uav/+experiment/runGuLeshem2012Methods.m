function record = runGuLeshem2012Methods(obs, cfg, methodNames)
%RUNGULESHEM2012METHODS 在同一论文例1观测上运行允许的比较方法。
%   允许方法仅包括原论文图1中出现且当前工程已经实现的Oracle、ESP、
%   WCO，以及本分支复现的gu_leshem2012_rec。所有非Oracle方法共享同一
%   SCM和名义期望导向矢量，真实量只在Oracle与指标评价阶段读取。

methodNames = reshape(string(methodNames), 1, []);
validateInputs(obs, cfg, methodNames);
paper = cfg.literature.guLeshem2012;
% 论文复现使用M=10；UPA迁移诊断复用同一比较器时从实际观测读取M。
M = size(obs.Rx, 1);
methodCfg = cfg;
methodCfg.derived.M = M;
methodCfg.beamforming.classic.esp.signalSubspaceDimension = ...
    paper.espSignalSubspaceDimension;
% 原论文WCO使用||e||^2<=0.3M；当前求解器的epsilonS表示||e||上界。
methodCfg.beamforming.epsilonS = sqrt(paper.wcoSquaredRadiusScale * M);

methodCount = numel(methodNames);
sinrLinear = zeros(1, methodCount);
desiredGain = zeros(1, methodCount);
interferenceLeakage = zeros(1, methodCount);
runtimeSec = zeros(1, methodCount);
incmNmse = nan(1, methodCount);
weights = complex(zeros(M, methodCount));
diagnostics = cell(1, methodCount);

for methodIndex = 1:methodCount
    timerValue = tic;
    method = methodNames(methodIndex);
    estimatedIncm = [];
    switch method
        case "oracle"
            [weight, beamformingDiagnostics] = uav.beamforming.solve( ...
                "mvdr", obs.truth.trueIncm, obs.truth.trueASVs(:, 1), methodCfg);
            methodDiagnostics = struct( ...
                'truthAccess', "oracle-only", ...
                'beamforming', beamformingDiagnostics);
        case "esp"
            [weight, beamformingDiagnostics] = ...
                uav.beamforming.eigenspaceProjection( ...
                obs.Rx, obs.nominalSteering, methodCfg);
            methodDiagnostics = struct('beamforming', beamformingDiagnostics);
        case "wco"
            [weight, beamformingDiagnostics] = uav.beamforming.solve( ...
                "wco", obs.Rx, obs.nominalSteering, methodCfg);
            methodDiagnostics = struct('beamforming', beamformingDiagnostics, ...
                'paperSquaredRadius', paper.wcoSquaredRadiusScale * M);
        case "gu_leshem2012_rec"
            [estimatedIncm, reconstructionDiagnostics] = ...
                uav.incm.guLeshem2012Rec( ...
                obs.Rx, obs.scanManifold, obs.scanWeights);
            [weight, beamformingDiagnostics] = uav.beamforming.solve( ...
                "mvdr", estimatedIncm, obs.nominalSteering, methodCfg);
            methodDiagnostics = struct( ...
                'reconstruction', reconstructionDiagnostics, ...
                'estimatedIncm', estimatedIncm, ...
                'beamforming', beamformingDiagnostics);
    end
    runtimeSec(methodIndex) = toc(timerValue);
    metricTruth = obs.truth;
    if ~isempty(estimatedIncm)
        metricTruth.estimatedIncm = estimatedIncm;
    end
    metrics = uav.metrics.evaluate(weight, metricTruth, methodCfg);
    sinrLinear(methodIndex) = metrics.sinrLinear;
    desiredGain(methodIndex) = metrics.desiredGain;
    interferenceLeakage(methodIndex) = metrics.interferenceLeakage;
    if isfield(metrics, 'incmNmse')
        incmNmse(methodIndex) = metrics.incmNmse;
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
    'weights', weights, ...
    'diagnostics', {diagnostics}, ...
    'snapshotSignatures', repmat(signature, methodCount, 1));
end

function validateInputs(obs, cfg, methodNames)
allowed = ["oracle", "esp", "wco", "gu_leshem2012_rec"];
if isempty(methodNames) || any(strlength(methodNames) == 0) || ...
        numel(unique(methodNames)) ~= numel(methodNames) || ...
        any(~ismember(methodNames, allowed))
    error('uav:experiment:InvalidGuLeshemMethodSelection', ...
        'Gu-Leshem comparison methods must be a unique nonempty allowed subset.');
end
if ~isstruct(cfg) || ~isfield(cfg, 'literature') || ...
        ~isfield(cfg.literature, 'guLeshem2012')
    error('uav:experiment:MissingGuLeshemConfiguration', ...
        'cfg.literature.guLeshem2012 is required.');
end
required = {'X', 'Rx', 'nominalSteering', 'scanManifold', ...
    'scanWeights', 'truth'};
if ~isstruct(obs) || ~all(isfield(obs, required))
    error('uav:experiment:InvalidGuLeshemObservation', ...
        'The paper-method runner requires a complete trial observation.');
end
end
