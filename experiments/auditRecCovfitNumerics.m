function audit = auditRecCovfitNumerics(resultFile, trialIndices, methodName)
%AUDITRECCOVFITNUMERICS 核查细化候选的拟合KKT及结果可复现性。
%   25°版本用已保存候选与系数独立重建；7°版本须重新生成候选并求解。
%   两种审计均不更换优化器或参数，不使用真值生成候选。
%   可读取运行中的检查点，但只审计其completedTrials，不把空白后缀算作结果。
arguments
    resultFile (1,1) string
    trialIndices (1,:) double = []
    methodName (1,1) string = "refined_candidate_covfit"
end
saved = load(resultFile);
if isfield(saved, 'results')
    r = saved.results;
    completed = size(r.sinrLinear, 1);
    isComplete = true;
else
    r = saved.checkpoint;
    completed = r.completedTrials;
    isComplete = false;
end
if isempty(trialIndices), trialIndices = 1:completed; end
assert(all(trialIndices >= 1 & trialIndices <= completed & mod(trialIndices, 1)==0), ...
    '只能审计已经完成的整数试次。');
cfg = r.cfg;
study = cfg.experiments.recRootCause;
assert(ismember(methodName,["refined_candidate_covfit","refined_narrow_sector"]), ...
    '数值审计只接受两种已记录的细化候选。');
methodIndex = find(r.methodNames == methodName);
rows = struct([]);
for trial = trialIndices
    for point = 1:numel(study.validationJitterDeg)
        local = cfg;
        local.attitude.stdDeg(:) = study.validationJitterDeg(point);
        scenario = uav.scenario.generate(local, uav.util.makeStream( ...
            cfg.random.baseSeed, study.validationSeedOffset+trial));
        for k = 1:study.validationPeriods
            obs = uav.scenario.generatePeriod(local, scenario, k, uav.util.makeStream( ...
                cfg.random.baseSeed, study.validationSeedOffset + ...
                study.validationTrials + (trial-1)*study.validationPeriods + k));
            for snapshotIndex = 1:numel(study.validationSnapshotCounts)
                L = study.validationSnapshotCounts(snapshotIndex);
                Rx = obs.X(:, 1:L) * obs.X(:, 1:L)' / L;
                if methodName == "refined_candidate_covfit"
                    d = r.refinementDiagnostics{trial, k, point, snapshotIndex};
                else
                    % 7°版本未逐点保存候选字典，须重放算法而不是借用25°诊断。
                    small = obs;
                    small.X = obs.X(:,1:L); small.Rx = Rx; small.pilot = obs.pilot(:,1:L);
                    candidate = r.candidateCfg;
                    candidate.timeline.L = L;
                    candidate.attitude.stdDeg = local.attitude.stdDeg;
                    candidate.spectrum.desiredExclusionHalfWidthDeg = study.narrowSectorHalfWidthDeg;
                    candidate.literature.guLeshem2012.upaDiagnosis.refineCandidates = true;
                    [~,d] = uav.incm.upaCovfit(uav.experiment.algorithmObservation(small),candidate);
                end
                assert(d.candidateRefinement.enabled && ~d.usedTruth, ...
                    '诊断数据不是启用细化且不读真值的方法。');
                A = d.candidateManifold;
                powers = d.fittedPowers;
                lambda = study.candidateCovfit.fitRegularization;
                residual = (A .* powers.') * A' + d.noisePower*eye(size(Rx)) - Rx;
                gradient = 2*real(sum(conj(A).*(residual*A), 1)).' + 2*lambda*powers;
                target = Rx - d.noisePower*eye(size(Rx));
                scale = max(norm(2*real(sum(conj(A).*(target*A), 1)), inf), eps);
                active = powers > sqrt(eps) * max(norm(powers, inf), 1);
                violations = [max(-powers, 0); max(-gradient, 0)/scale; ...
                    abs(gradient(active))/scale];
                desiredCount = size(d.desiredCandidateNodes, 1);
                Ai = A(:, desiredCount+1:end);
                pi = powers(desiredCount+1:end);
                Ri = (Ai .* pi.') * Ai' + d.noisePower*eye(size(Rx));
                w = uav.beamforming.solve("mvdr", Ri, cfg.derived.nominalSteeringReference, cfg);
                m = uav.metrics.evaluate(w, obs.truth, cfg);
                reference = r.sinrLinear(trial, k, point, snapshotIndex, methodIndex);
                flags = [d.candidateRefinement.desired.exitflags; ...
                    d.candidateRefinement.interferer.exitflags];
                accepted = [d.candidateRefinement.desired.accepted; ...
                    d.candidateRefinement.interferer.accepted];
                row = struct('trial', trial, 'period', k, 'jitterDeg', study.validationJitterDeg(point), ...
                    'L', L, 'kktRelativeViolation', max(violations), ...
                    'relativeSinrReplayError', abs(m.sinrLinear-reference)/max(reference, eps), ...
                    'refinementAttemptCount', numel(flags), ...
                    'nonpositiveExitCount', sum(flags<=0), 'fallbackCount', sum(~accepted), ...
                    'candidateCount', numel(powers), 'fitResidual', d.fitResidual);
                rows = [rows; row]; %#ok<AGROW>
            end
        end
    end
    fprintf('NUMERICAL_AUDIT_PROGRESS target=%s trial=%d/%d\n',methodName,trial,completed);
end
audit = struct('sourceFile', resultFile, 'sourceIsComplete', isComplete, ...
    'methodName',methodName, ...
    'completedTrialsAtLoad', completed, 'auditedTrials', trialIndices, ...
    'records', struct2table(rows), 'generatedAt', datetime('now'));
directory = fileparts(resultFile);
tag = "numerical_audit_" + cfg.profile;
if methodName == "refined_narrow_sector", tag = tag+"_"+methodName; end
if ~isComplete, tag = tag + "_partial"; end
save(fullfile(directory, tag + ".mat"), 'audit');
writetable(audit.records, fullfile(directory, tag + ".csv"));
fprintf('NUMERICAL_AUDIT trials=%d points=%d maxKKT=%.6g maxReplay=%.6g fallbacks=%d/%d\n', ...
    numel(trialIndices), height(audit.records), max(audit.records.kktRelativeViolation), ...
    max(audit.records.relativeSinrReplayError), sum(audit.records.fallbackCount), ...
    sum(audit.records.refinementAttemptCount));
end
