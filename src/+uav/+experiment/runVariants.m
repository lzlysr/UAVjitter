function [record, nextState] = runVariants(obs, priorState, cfg)
%RUNVARIANTS Evaluate full temporal reconstruction and one-switch ablations.

variantNames = reshape(cfg.experiments.ablationNames, 1, []);
expected = ["full", "lambda_t_zero", "identity_alignment", ...
    "uniform_quadrature", "no_structure", "pilot_projection"];
if ~isequal(variantNames, expected)
    error('uav:experiment:InvalidAblationContract', ...
        'The ablation runner requires the fixed configured variant order.');
end
if isempty(priorState)
    priorState = cell(1, numel(variantNames));
elseif ~iscell(priorState) || numel(priorState) ~= numel(variantNames)
    error('uav:experiment:InvalidAblationState', ...
        'Ablation state must contain one independent temporal state per variant.');
end

switchDeltaFields = ["none", "incm.lambdaT", "incm.ablation.alignment", ...
    "incm.ablation.adaptiveQuadrature", ...
    "incm.ablation.structureEnhancement", ...
    "incm.ablation.pilotProjection"];
variantCount = numel(variantNames);
sinrLinear = zeros(1, variantCount);
desiredGain = zeros(1, variantCount);
interferenceLeakage = zeros(1, variantCount);
runtimeSec = zeros(1, variantCount);
incmNmse = zeros(1, variantCount);
diagnostics = cell(1, variantCount);
effectiveConfigs = cell(1, variantCount);
aRef = cfg.derived.nominalSteeringReference;
for variantIndex = 1:variantCount
    variantCfg = applyVariant(cfg, variantNames(variantIndex));
    effectiveConfigs{variantIndex} = variantCfg.incm;
    variantTimer = tic;
    [estimatedIncm, priorState{variantIndex}, reconstructionDiagnostics] = ...
        uav.incm.reconstruct("temporal", obs, priorState{variantIndex}, variantCfg);
    [weight, beamformingDiagnostics] = uav.beamforming.solve( ...
        "robust", estimatedIncm, aRef, variantCfg);
    runtimeSec(variantIndex) = toc(variantTimer);
    metricTruth = obs.truth;
    metricTruth.estimatedIncm = estimatedIncm;
    metrics = uav.metrics.evaluate(weight, metricTruth, variantCfg);
    sinrLinear(variantIndex) = metrics.sinrLinear;
    desiredGain(variantIndex) = metrics.desiredGain;
    interferenceLeakage(variantIndex) = metrics.interferenceLeakage;
    incmNmse(variantIndex) = metrics.incmNmse;
    diagnostics{variantIndex} = struct( ...
        'reconstruction', reconstructionDiagnostics, ...
        'beamforming', beamformingDiagnostics, ...
        'changedField', switchDeltaFields(variantIndex));
end
signature = uav.util.snapshotSignature(obs.X);
record = struct( ...
    'variantNames', variantNames, ...
    'switchDeltaFields', switchDeltaFields, ...
    'sinrLinear', sinrLinear, ...
    'desiredGain', desiredGain, ...
    'interferenceLeakage', interferenceLeakage, ...
    'incmNmse', incmNmse, ...
    'runtimeSec', runtimeSec, ...
    'diagnostics', {diagnostics}, ...
    'effectiveIncmConfigs', {effectiveConfigs}, ...
    'snapshotSignatures', repmat(signature, variantCount, 1), ...
    'epsilonS', cfg.beamforming.epsilonS, ...
    'gamma', cfg.beamforming.gamma);
nextState = priorState;
end

function cfg = applyVariant(cfg, variant)
switch variant
    case "full"
    case "lambda_t_zero"
        cfg.incm.lambdaT = 0;
    case "identity_alignment"
        cfg.incm.ablation.alignment = false;
    case "uniform_quadrature"
        cfg.incm.ablation.adaptiveQuadrature = false;
    case "no_structure"
        cfg.incm.ablation.structureEnhancement = false;
    case "pilot_projection"
        cfg.incm.ablation.pilotProjection = true;
end
end
