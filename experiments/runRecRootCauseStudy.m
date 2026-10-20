function results = runRecRootCauseStudy(profile)
%RUNRECROOTCAUSESTUDY 在零失配下隔离原始Rec的样本、测度及背景影响。
%   精确总协方差仅是诊断干预，不是向可部署算法泄露真值。
%   raw始终调用原论文公式核心；其他变体仅用于因果辨别，不能冒充原式。
arguments
    profile (1,1) string = "smoke"
end
cfg = config(profile);
cfg.attitude.stdDeg = zeros(size(cfg.attitude.stdDeg));
cfg = uav.config.validate(cfg);
study = cfg.experiments.recRootCause;
studyTimer = tic;
stream = RandStream(char(cfg.random.generator), 'Seed', cfg.random.baseSeed);
scenario = uav.scenario.generate(cfg, stream);
obs = uav.scenario.generatePeriod(cfg, scenario, 1, stream);
A = obs.truth.trueASVs;
Rin = obs.truth.trueIncm;
R = Rin + cfg.derived.desiredPower * A(:, 1) * A(:, 1)';
records = struct([]);
% 网格收敛：其余模型、角域、输入矩阵和权值接口完全一致。
for measure = study.measures
    for gridSize = study.gridSizes
        grid = uav.spectrum.upaIntegrationGrid(cfg, gridSize, measure);
        records = [records; examine(R, Rin, A(:, 1), grid, cfg, ...
            "population-grid", gridSize, "raw", 0)]; %#ok<AGROW>
    end
end
% 角域和谱分解干预：纯噪声Capon值sigma/M只作为解析诊断参照。
for halfWidth = study.exclusionHalfWidthsDeg
    local = cfg;
    local.spectrum.desiredExclusionHalfWidthDeg = halfWidth;
    grid = uav.spectrum.upaIntegrationGrid(local, study.diagnosticGridSize, "solid_angle");
    records = [records; examine(R, Rin, A(:, 1), grid, local, ...
        "exclusion", halfWidth, "raw", 0)]; %#ok<AGROW>
end
grid = uav.spectrum.upaIntegrationGrid(cfg, study.diagnosticGridSize, "solid_angle");
for variant = ["raw", "remove-floor", "remove-floor-add-white", ...
        "floor-only", "remove-soi-from-input", "true-interference-white", ...
        "known-directions-capon-height", "raw-model-projection", ...
        "whiten-floor-trace-preserved", "collapse-excess-true-directions", ...
        "collapse-and-whiten-floor"]
    records = [records; examine(R, Rin, A(:, 1), grid, cfg, ...
        "components", 0, variant, 0)]; %#ok<AGROW>
end
% 加载、谱底、白噪声项的全因子干预：全部共享精确输入及同一网格。
for loadingRatio = [0, cfg.spectrum.diagonalLoading]
    for subtractFloor = [false, true]
        for addWhite = [false, true]
            intervention = struct('loadingRatio', loadingRatio, ...
                'subtractFloor', subtractFloor, 'addWhite', addWhite);
            label = sprintf('load=%g-floor=%d-white=%d', ...
                loadingRatio, subtractFloor, addWhite);
            records = [records; examine(R, Rin, A(:, 1), grid, cfg, ...
                "engineering-factorial", loadingRatio, string(label), 0, intervention)]; %#ok<AGROW>
        end
    end
end
% 将原始矩阵精确分为已知点源模型投影与残差，只缩放其中一个部分。
% 这些干预读取真方向，仅用于定位损失，绝不是可部署方法。
for scale = study.matrixResidualScales
    records = [records; examine(R, Rin, A(:, 1), grid, cfg, ...
        "matrix-residual-scale", scale, "diagnostic-truth-dictionary", 0)]; %#ok<AGROW>
end
for scale = study.pointPowerScales
    records = [records; examine(R, Rin, A(:, 1), grid, cfg, ...
        "point-power-scale", scale, "diagnostic-truth-dictionary", 0)]; %#ok<AGROW>
end
% 只改变干扰功率；不改变参考方向或利用真值选择积分区域。
for scale = study.interferencePowerScales
    Ri = cfg.derived.noisePower * eye(cfg.derived.M) + ...
        scale * (Rin - cfg.derived.noisePower * eye(cfg.derived.M));
    Rx = Ri + cfg.derived.desiredPower * A(:, 1) * A(:, 1)';
    records = [records; examine(Rx, Ri, A(:, 1), grid, cfg, ...
        "interference-scale", scale, "raw", 0)]; %#ok<AGROW>
end
% 同为二维UPA，分别改变孔径或源几何，检查“二维结构必然失效”的断言。
for shapeIndex = 1:size(study.arrayShapes, 1)
    local = cfg;
    local.array.Nx = study.arrayShapes(shapeIndex, 1);
    local.array.Ny = study.arrayShapes(shapeIndex, 2);
    local = uav.config.validate(local);
    [exactRx, exactRi, exactA] = exactNominalCovariance(local);
    localGrid = uav.spectrum.upaIntegrationGrid(local, ...
        study.diagnosticGridSize, "solid_angle");
    records = [records; examine(exactRx, exactRi, exactA(:, 1), localGrid, ...
        local, "array-aperture", local.derived.M, "raw", 0)]; %#ok<AGROW>
end
local = cfg;
local.sources.desiredAzElDeg = study.orthogonalDesiredAzElDeg;
local.sources.interfererAzElDeg = study.orthogonalInterfererAzElDeg;
local = uav.config.validate(local);
[exactRx, exactRi, exactA] = exactNominalCovariance(local);
for gridSize = study.gridSizes
    localGrid = uav.spectrum.upaIntegrationGrid(local, gridSize, "solid_angle");
    records = [records; examine(exactRx, exactRi, exactA(:, 1), localGrid, ...
        local, "orthogonal-geometry", gridSize, "raw", 0)]; %#ok<AGROW>
end
% 同一批最大长度快拍取前缀，各方法及各L严格共享观测实现。
for trial = 1:study.trialCount
    stream = RandStream(char(cfg.random.generator), 'Seed', cfg.random.baseSeed);
    stream.Substream = trial;
    local = cfg;
    local.timeline.L = max(study.snapshotCounts);
    observed = uav.scenario.generatePeriod(local, scenario, 1, stream);
    for L = study.snapshotCounts
        Rx = observed.X(:, 1:L) * observed.X(:, 1:L)' / L;
        if L < cfg.derived.M
            % 原式要求可逆，不能用加载或伪逆偷偷补齐。
            record = emptyRecord("snapshots", L, "raw", trial);
            record.status = "not-applicable-L<M";
        else
            record = examine(Rx, Rin, A(:, 1), grid, cfg, ...
                "snapshots", L, "raw", trial);
        end
        records = [records; record]; %#ok<AGROW>
    end
end
root = fileparts(fileparts(mfilename('fullpath')));
outputDirectory = fullfile(root, study.outputDirectory);
if ~isfolder(outputDirectory), mkdir(outputDirectory); end
results = struct('cfg', cfg, 'records', struct2table(records), ...
    'populationCovariance', R, 'trueIncm', Rin, 'trueManifold', A, ...
    'metadata', uav.experiment.metadata(cfg, "rec-root-cause", ...
        struct('trials', study.trialCount), toc(studyTimer)), ...
    'experimentSource', fileread(mfilename('fullpath') + ".m"), ...
    'configurationSource', fileread(fullfile(root, 'config.m')), ...
    'generatedAt', datetime('now'), 'diagnosticOnlyTruth', true);
save(fullfile(outputDirectory, 'root_cause.mat'), 'results', '-v7.3');
writetable(results.records, fullfile(outputDirectory, 'root_cause.csv'));
disp(results.records(results.records.stage ~= "snapshots", :));
fprintf('ROOT_CAUSE_OUTPUT %s\n', outputDirectory);
end

function [Rx, Ri, A] = exactNominalCovariance(cfg)
%EXACTNOMINALCOVARIANCE 仅供诊断的解析输入，不生成额外样本误差。
A = cfg.derived.nominalSteering;
Ri = cfg.derived.noisePower * eye(cfg.derived.M) + ...
    A(:, 2:end) * diag(cfg.derived.interfererPowers) * A(:, 2:end)';
Rx = Ri + cfg.derived.desiredPower * A(:, 1) * A(:, 1)';
end

function record = examine(Rx, truth, a, grid, cfg, stage, value, variant, trial, intervention)
%EXAMINE 记录SINR、白噪声放大及干扰泄漏，避免仅看一个最终指标。
record = emptyRecord(stage, value, variant, trial);
record.measure = grid.measureName;
mask = grid.includedMask;
G = grid.manifold(:, mask);
weights = grid.weights(mask);
if variant == "remove-soi-from-input", Rx = truth; end
[Rhat, diagnostics] = uav.incm.guLeshem2012Rec(Rx, G, weights);
floorPower = cfg.derived.noisePower / cfg.derived.M;
floorCovariance = (G .* (weights * floorPower)) * G';
if nargin >= 10
    loading = intervention.loadingRatio * real(trace(Rx)) / size(Rx, 1);
    [Rhat, diagnostics] = uav.incm.guLeshem2012Rec( ...
        Rx + loading * eye(size(Rx)), G, weights);
    loadedFloor = (G .* (weights * (cfg.derived.noisePower + loading) / ...
        cfg.derived.M)) * G';
    if intervention.subtractFloor, Rhat = Rhat - loadedFloor; end
    if intervention.addWhite, Rhat = Rhat + cfg.derived.noisePower * eye(size(Rhat)); end
end
switch variant
    case "remove-floor"
        Rhat = Rhat - floorCovariance;
    case "remove-floor-add-white"
        Rhat = Rhat - floorCovariance + cfg.derived.noisePower * eye(size(Rhat));
    case "floor-only"
        Rhat = floorCovariance;
    case "true-interference-white"
        Rhat = truth;
end
Rhat = (Rhat + Rhat') / 2;
% 已知方向仅用于诊断：检验面积积分后的矩阵能否被点源模型解释。
Ai = cfg.derived.nominalSteering(:, 2:end);
dictionary = zeros(numel(Rhat), size(Ai, 2) + 1);
for source = 1:size(Ai, 2)
    covariance = Ai(:, source) * Ai(:, source)';
    dictionary(:, source) = covariance(:);
end
identity = eye(size(Rhat));
dictionary(:, end) = identity(:);
modelPowers = lsqnonneg([real(dictionary); imag(dictionary)], ...
    [real(Rhat(:)); imag(Rhat(:))]);
modelCovariance = reshape(dictionary * modelPowers, size(Rhat));
record.pointModelResidual = norm(Rhat - modelCovariance, 'fro') / norm(Rhat, 'fro');
record.fittedNoise = modelPowers(end);
record.fittedInterferencePowerSum = sum(modelPowers(1:end-1));
if variant == "raw-model-projection"
    Rhat = modelCovariance;
elseif variant == "known-directions-capon-height"
    caponHeight = 1 ./ real(sum(conj(Ai) .* (Rx \ Ai), 1));
    pointPowers = max(caponHeight - floorPower, 0);
    Rhat = (Ai .* pointPowers) * Ai' + cfg.derived.noisePower * eye(size(Rhat));
elseif variant == "whiten-floor-trace-preserved"
    Rhat = Rhat - floorCovariance + real(trace(floorCovariance)) / size(Rhat, 1) * eye(size(Rhat));
elseif variant == "collapse-excess-true-directions" || variant == "collapse-and-whiten-floor"
    % 保持谱积分的总质量不变，只将各谱节点指纹替换为最近真实干扰指纹。
    % 此项读真方向，仅作因果干预，绝非候选筛选或部署算法。
    [~, closest] = max(abs(Ai' * G).^2, [], 1);
    mass = accumarray(closest.', ...
        (weights .* max(diagnostics.caponPower-floorPower, 0)).', [size(Ai, 2), 1]);
    Rhat = Ai * diag(mass) * Ai' + floorCovariance;
    if variant == "collapse-and-whiten-floor"
        Rhat = Ai * diag(mass) * Ai' + ...
            real(trace(floorCovariance)) / size(Rhat, 1) * eye(size(Rhat));
    end
end
if stage == "matrix-residual-scale"
    Rhat = modelCovariance + value * (Rhat - modelCovariance);
elseif stage == "point-power-scale"
    Rhat = value * Ai * diag(modelPowers(1:end-1)) * Ai' + ...
        modelPowers(end) * eye(size(Rhat));
end
[w, wd] = uav.beamforming.solve("mvdr", Rhat, a, cfg);
[wo, ~] = uav.beamforming.solve("mvdr", truth, a, cfg);
signalPower = cfg.derived.desiredPower * abs(w' * a)^2;
record.sinrDb = 10 * log10(signalPower / real(w' * truth * w));
record.oracleSinrDb = 10 * log10(cfg.derived.desiredPower * ...
    abs(wo' * a)^2 / real(wo' * truth * wo));
record.noiseOutput = cfg.derived.noisePower * real(w' * w);
record.interferenceOutput = real(w' * (truth - ...
    cfg.derived.noisePower * eye(size(truth))) * w);
record.traceRatio = real(trace(Rhat) / trace(truth));
record.floorTraceFraction = real(trace(floorCovariance)) / ...
    sum(diagnostics.weightedPower) / cfg.derived.M;
record.minimumEigenvalue = min(real(eig(Rhat)));
record.weightLoading = wd.loading;
record.status = "ok";
if trial == 0
    fprintf('ROOT_CAUSE %s %.3g %s %s SINR=%.6f noise=%.6g interference=%.6g\n', ...
        stage, value, grid.measureName, variant, record.sinrDb, ...
        record.noiseOutput, record.interferenceOutput);
end
end

function record = emptyRecord(stage, value, variant, trial)
%EMPTYRECORD 对不适用点保留显式状态和NaN，不伪造原式数值。
record = struct('stage', string(stage), 'value', value, 'variant', string(variant), ...
    'trial', trial, 'measure', "", 'sinrDb', NaN, 'oracleSinrDb', NaN, ...
    'noiseOutput', NaN, 'interferenceOutput', NaN, 'traceRatio', NaN, ...
    'floorTraceFraction', NaN, 'minimumEigenvalue', NaN, 'weightLoading', NaN, ...
    'pointModelResidual', NaN, 'fittedNoise', NaN, ...
    'fittedInterferencePowerSum', NaN, ...
    'status', "pending");
end
