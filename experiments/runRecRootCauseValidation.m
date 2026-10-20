function results = runRecRootCauseValidation(profile)
%RUNRECROOTCAUSEVALIDATION 用独立配对子流验证点源模型候选，不改变鲁棒基线。
%   candidate_covfit仅是已有CovFit的预先固定参数变体，不是原始Rec。
%   Oracle和true_incm_nominal均为真值诊断，不能算可部署方法。
arguments
    profile (1,1) string = "smoke"
end
cfg = config(profile);
study = cfg.experiments.recRootCause;
% 保留配置中其他实验的合法K上界，本实验只消费前validationPeriods周期。
cfg.timeline.K = max(cfg.timeline.K, study.validationPeriods);
cfg.timeline.L = max(study.validationSnapshotCounts);
cfg.spectrum.gridSize = cfg.literature.guLeshem2012.upaDiagnosis.selectedGridSize;
cfg = uav.config.validate(cfg);
candidate = cfg;
candidate.spectrum.gridSize = study.candidateCovfit.gridSize;
candidate.literature.guLeshem2012.upaDiagnosis.loadingRatio = study.candidateCovfit.loadingRatio;
candidate.literature.guLeshem2012.upaDiagnosis.peakNmsRadiusBins = study.candidateCovfit.peakNmsRadiusBins;
candidate.literature.guLeshem2012.upaDiagnosis.fitRegularization = study.candidateCovfit.fitRegularization;
candidate = uav.config.validate(candidate);
grid = uav.spectrum.upaIntegrationGrid(cfg, study.diagnosticGridSize, "solid_angle");
baseNames = ["oracle", "esp", "wco", "gu_leshem2012_rec", "upa_covfit_incm"];
methodNames = [baseNames, "candidate_covfit", "dl_mvdr", "old_static", ...
    "true_incm_nominal", "refined_candidate_covfit", "refined_narrow_sector"];
dimensions = [study.validationTrials, study.validationPeriods, ...
    numel(study.validationJitterDeg), numel(study.validationSnapshotCounts), numel(methodNames)];
sinrLinear = nan(dimensions);
desiredGain = nan(dimensions);
interferenceLeakage = nan(dimensions);
incmNmse = nan(dimensions);
runtimeSec = nan(dimensions);
root = fileparts(fileparts(mfilename('fullpath')));
directory = fullfile(root, study.outputDirectory);
if ~isfolder(directory), mkdir(directory); end
refinementDiagnostics = cell(dimensions(1:4));
desiredOutsidePrior = false(dimensions(1:3));
interfererInsidePrior = false(dimensions(1:3));
timer = tic;
for trial = 1:study.validationTrials
    for jitterIndex = 1:numel(study.validationJitterDeg)
        local = cfg;
        local.attitude.stdDeg(:) = study.validationJitterDeg(jitterIndex);
        % 各抖动点共用标准正态实现，改变的只是RMS缩放。
        scenario = uav.scenario.generate(local, uav.util.makeStream( ...
            cfg.random.baseSeed, study.validationSeedOffset + trial));
        for k = 1:study.validationPeriods
            obs = uav.scenario.generatePeriod(local, scenario, k, uav.util.makeStream( ...
                cfg.random.baseSeed, study.validationSeedOffset + ...
                study.validationTrials + (trial-1)*study.validationPeriods + k));
            % 只用于事后归因：复现算法下半球节点映射后的真实源分区。
            truthOutside = uav.spectrum.admissibleDirections(obs.truth.directionCosines.', local);
            desiredOutsidePrior(trial, k, jitterIndex) = truthOutside(1);
            interfererInsidePrior(trial, k, jitterIndex) = any(~truthOutside(2:end));
            for snapshotIndex = 1:numel(study.validationSnapshotCounts)
                L = study.validationSnapshotCounts(snapshotIndex);
                small = obs;
                small.X = obs.X(:, 1:L);
                small.Rx = small.X * small.X' / L;
                small.pilot = obs.pilot(:, 1:L);
                pointCfg = local;
                pointCfg.timeline.L = L;
                pointCandidate = candidate;
                pointCandidate.timeline.L = L;
                pointCandidate.attitude.stdDeg = local.attitude.stdDeg;
                record = uav.experiment.runUpaDiagnosisMethods(small, pointCfg, baseNames, grid);
                values = record.sinrLinear;
                gains = record.desiredGain;
                leaks = record.interferenceLeakage;
                errors = record.incmNmse;
                times = record.runtimeSec;
                safe = uav.experiment.algorithmObservation(small);
                extraTimer = tic;
                Rh = uav.incm.upaCovfit(safe, pointCandidate);
                [m, elapsed] = evaluateMatrix(Rh, safe.nominalSteeringReference, ...
                    small.truth, pointCfg, "mvdr", extraTimer);
                values(6) = m.sinrLinear; gains(6) = m.desiredGain;
                leaks(6) = m.interferenceLeakage; errors(6) = m.incmNmse; times(6) = elapsed;
                extraTimer = tic;
                [m, elapsed] = evaluateMatrix(safe.Rx, safe.nominalSteeringReference, ...
                    small.truth, pointCfg, "dlmvdr", extraTimer);
                values(7) = m.sinrLinear; gains(7) = m.desiredGain;
                leaks(7) = m.interferenceLeakage; errors(7) = NaN; times(7) = elapsed;
                extraTimer = tic;
                Rh = uav.incm.reconstruct("static", safe, [], pointCfg);
                [m, elapsed] = evaluateMatrix(Rh, safe.nominalSteeringReference, ...
                    small.truth, pointCfg, "mvdr", extraTimer);
                values(8) = m.sinrLinear; gains(8) = m.desiredGain;
                leaks(8) = m.interferenceLeakage; errors(8) = m.incmNmse; times(8) = elapsed;
                extraTimer = tic;
                [m, elapsed] = evaluateMatrix(small.truth.trueIncm, safe.nominalSteeringReference, ...
                    small.truth, pointCfg, "mvdr", extraTimer);
                values(9) = m.sinrLinear; gains(9) = m.desiredGain;
                leaks(9) = m.interferenceLeakage; errors(9) = 0; times(9) = elapsed;
                extraTimer = tic;
                pointCandidate.literature.guLeshem2012.upaDiagnosis.refineCandidates = true;
                [Rh, refinedDiagnostic] = uav.incm.upaCovfit(safe, pointCandidate);
                refinementDiagnostics{trial, k, jitterIndex, snapshotIndex} = refinedDiagnostic;
                [m, elapsed] = evaluateMatrix(Rh, safe.nominalSteeringReference, ...
                    small.truth, pointCfg, "mvdr", extraTimer);
                values(10) = m.sinrLinear; gains(10) = m.desiredGain;
                leaks(10) = m.interferenceLeakage; errors(10) = m.incmNmse; times(10) = elapsed;
                extraTimer = tic;
                pointCandidate.spectrum.desiredExclusionHalfWidthDeg = study.narrowSectorHalfWidthDeg;
                Rh = uav.incm.upaCovfit(safe, pointCandidate);
                [m, elapsed] = evaluateMatrix(Rh, safe.nominalSteeringReference, ...
                    small.truth, pointCfg, "mvdr", extraTimer);
                values(11) = m.sinrLinear; gains(11) = m.desiredGain;
                leaks(11) = m.interferenceLeakage; errors(11) = m.incmNmse; times(11) = elapsed;
                sinrLinear(trial, k, jitterIndex, snapshotIndex, :) = values;
                desiredGain(trial, k, jitterIndex, snapshotIndex, :) = gains;
                interferenceLeakage(trial, k, jitterIndex, snapshotIndex, :) = leaks;
                incmNmse(trial, k, jitterIndex, snapshotIndex, :) = errors;
                runtimeSec(trial, k, jitterIndex, snapshotIndex, :) = times;
            end
        end
    end
    fprintf('PAIRED_VALIDATION trial=%d/%d elapsed=%.1fs\n', trial, study.validationTrials, toc(timer));
    % 每个独立试次落盘，长验证即使中断也能保留已完成原始数据。
    checkpoint = struct('cfg', cfg, 'candidateCfg', candidate, 'methodNames', methodNames, ...
        'completedTrials', trial, 'sinrLinear', sinrLinear, 'desiredGain', desiredGain, ...
        'interferenceLeakage', interferenceLeakage, 'incmNmse', incmNmse, ...
        'runtimeSec', runtimeSec, 'refinementDiagnostics', {refinementDiagnostics}, ...
        'desiredOutsidePrior', desiredOutsidePrior, 'interfererInsidePrior', interfererInsidePrior);
    save(fullfile(directory, "paired_validation_" + profile + "_checkpoint.mat"), 'checkpoint', '-v7.3');
end
summaryRecords = struct([]);
for snapshotIndex = 1:numel(study.validationSnapshotCounts)
    for jitterIndex = 1:numel(study.validationJitterDeg)
        for methodIndex = 1:numel(methodNames)
            raw = sinrLinear(:, :, jitterIndex, snapshotIndex, methodIndex);
            % 先对每个独立试次在线性域平均周期，再在试次间取分位。
            trialMean = mean(raw, 2);
            summary = struct('L', study.validationSnapshotCounts(snapshotIndex), ...
                'jitterDeg', study.validationJitterDeg(jitterIndex), ...
                'method', methodNames(methodIndex), ...
                'medianSinrDb', 10*log10(median(trialMean)), ...
                'lowSinrDb', 10*log10(prctile(trialMean, 5)), ...
                'periodLowSinrDb', 10*log10(prctile(raw(:), 5)), ...
                'outageProbability', sum(raw(:) < 10^(study.validationOutageThresholdDb/10)) / sum(isfinite(raw(:))), ...
                'medianIncmNmse', median(incmNmse(:, :, jitterIndex, snapshotIndex, methodIndex), 'all'));
            summaryRecords = [summaryRecords; summary]; %#ok<AGROW>
        end
    end
end
% 以独立试次为重采样单位，同一抽样索引用于两种方法；不把相关周期当独立样本。
pairedRecords = struct([]);
comparisonNames = ["esp", "wco", "dl_mvdr", "upa_covfit_incm"];
comparisonStream = uav.util.makeStream(cfg.random.baseSeed, study.validationSeedOffset*2);
for snapshotIndex = 1:numel(study.validationSnapshotCounts)
    for jitterIndex = 1:numel(study.validationJitterDeg)
        for targetIndex = [10, 11]
            targetTrials = mean(sinrLinear(:, :, jitterIndex, snapshotIndex, targetIndex), 2);
            for comparison = comparisonNames
                referenceIndex = find(methodNames == comparison);
                referenceTrials = mean(sinrLinear(:, :, jitterIndex, snapshotIndex, referenceIndex), 2);
                indices = randi(comparisonStream, study.validationTrials, ...
                    study.validationTrials, cfg.experiments.bootstrapReplicates);
                differences = 10*log10(median(targetTrials(indices), 1)) - ...
                    10*log10(median(referenceTrials(indices), 1));
                interval = prctile(differences, 100*cfg.experiments.aggregation.confidenceQuantiles);
                pair = struct('L', study.validationSnapshotCounts(snapshotIndex), ...
                    'jitterDeg', study.validationJitterDeg(jitterIndex), ...
                    'target', methodNames(targetIndex), 'reference', comparison, ...
                    'medianDifferenceDb', 10*log10(median(targetTrials)/median(referenceTrials)), ...
                    'ciLowerDb', interval(1), 'ciUpperDb', interval(2));
                pairedRecords = [pairedRecords; pair]; %#ok<AGROW>
            end
        end
    end
end
results = struct('cfg', cfg, 'candidateCfg', candidate, 'methodNames', methodNames, ...
    'sinrLinear', sinrLinear, 'desiredGain', desiredGain, ...
    'interferenceLeakage', interferenceLeakage, 'incmNmse', incmNmse, ...
    'runtimeSec', runtimeSec, 'summary', struct2table(summaryRecords), ...
    'pairedComparisons', struct2table(pairedRecords), ...
    'refinementDiagnostics', {refinementDiagnostics}, ...
    'desiredOutsidePrior', desiredOutsidePrior, ...
    'interfererInsidePrior', interfererInsidePrior, ...
    'metadata', uav.experiment.metadata(cfg, "root-cause-paired-validation", ...
        struct('trials', study.validationTrials, 'periods', study.validationPeriods), toc(timer)), ...
    'experimentSource', fileread(mfilename('fullpath') + ".m"));
results.sourceSnapshot = struct('config', fileread(fullfile(root, 'config.m')), ...
    'covfit', fileread(fullfile(root, 'src', '+uav', '+incm', 'upaCovfit.m')), ...
    'refinement', fileread(fullfile(root, 'src', '+uav', '+spectrum', 'refineCaponNodes.m')));
save(fullfile(directory, "paired_validation_" + profile + ".mat"), 'results', '-v7.3');
writetable(results.summary, fullfile(directory, "paired_validation_" + profile + ".csv"));
writetable(results.pairedComparisons, fullfile(directory, "paired_comparisons_" + profile + ".csv"));
disp(results.summary);
end

function [metrics, elapsed] = evaluateMatrix(R, aRef, truth, cfg, method, timer)
%EVALUATEMATRIX 先完成权值求解并计时，再单独用真值评价，真值不进入求解器。
w = uav.beamforming.solve(method, R, aRef, cfg);
elapsed = toc(timer);
truth.estimatedIncm = R;
metrics = uav.metrics.evaluate(w, truth, cfg);
end
