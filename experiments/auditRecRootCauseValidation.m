function audit = auditRecRootCauseValidation(resultFile)
%AUDITRECROOTCAUSEVALIDATION 对完整验证作覆盖、真值上界及最差周期复核。
%   此函数不接受检查点作为完整证据，不重新选参数或删除失败点。
arguments
    resultFile (1,1) string
end
saved = load(resultFile, 'results');
assert(isfield(saved, 'results'), '必须提供完整结果，不能把检查点当作完成证据。');
r = saved.results;
cfg = r.cfg;
s = cfg.experiments.recRootCause;
root = fileparts(fileparts(mfilename('fullpath')));
assert(strcmp(r.sourceSnapshot.config,fileread(fullfile(root,'config.m'))) && ...
    strcmp(r.sourceSnapshot.covfit,fileread(fullfile(root,'src','+uav','+incm','upaCovfit.m'))) && ...
    strcmp(r.sourceSnapshot.refinement,fileread(fullfile(root,'src','+uav','+spectrum','refineCaponNodes.m'))) && ...
    strcmp(r.experimentSource,fileread(fullfile(root,'experiments','runRecRootCauseValidation.m'))), ...
    '当前配置或仿真算法源码与结果中的快照不一致，不能直接归为同一版本。');
assert(size(r.sinrLinear, 1)==s.validationTrials && ...
    size(r.sinrLinear, 2)==s.validationPeriods && ...
    size(r.sinrLinear, 3)==numel(s.validationJitterDeg) && ...
    size(r.sinrLinear, 4)==numel(s.validationSnapshotCounts) && ...
    size(r.sinrLinear, 5)==numel(r.methodNames), '实际维度与验证协议不一致。');
originalIndex = find(r.methodNames == "gu_leshem2012_rec");
allowedMissing = false(size(r.sinrLinear));
allowedMissing(:, :, :, s.validationSnapshotCounts < cfg.derived.M, originalIndex) = true;
assert(all(isnan(r.sinrLinear(allowedMissing))), 'L<M的原式不能有伪造成绩。');
assert(all(isfinite(r.sinrLinear(~allowedMissing)) & r.sinrLinear(~allowedMissing)>0), ...
    '适用点存在非有限或非正SINR。');
oracle = r.sinrLinear(:, :, :, :, r.methodNames=="oracle");
excess = r.sinrLinear - oracle;
tolerance = cfg.experiments.oracleTolerance * max(oracle, 1);
expandedTolerance = repmat(tolerance, 1,1,1,1,numel(r.methodNames));
assert(all(excess(~allowedMissing) <= expandedTolerance(~allowedMissing)), ...
    '存在超越真实方向Oracle的权值，须检查指标或数据流。');
% 分区统计使用重放的真方向，仅供结果解释，绝不反馈候选筛选。
priorRecords = struct([]);
for halfWidth = [cfg.spectrum.desiredExclusionHalfWidthDeg, s.narrowSectorHalfWidthDeg]
    for point = 1:numel(s.validationJitterDeg)
        desiredOutside = false(s.validationTrials, s.validationPeriods);
        interfererInside = desiredOutside;
        for trial = 1:s.validationTrials
            local = cfg;
            local.attitude.stdDeg(:) = s.validationJitterDeg(point);
            local.spectrum.desiredExclusionHalfWidthDeg = halfWidth;
            scene = uav.scenario.generate(local, uav.util.makeStream( ...
                cfg.random.baseSeed, s.validationSeedOffset+trial));
            for k = 1:s.validationPeriods
                outside = uav.spectrum.admissibleDirections(scene.directionCosines(:, :, k).', local);
                desiredOutside(trial, k) = outside(1);
                interfererInside(trial, k) = any(~outside(2:end));
            end
        end
        row = struct('halfWidthDeg', halfWidth, 'jitterDeg', s.validationJitterDeg(point), ...
            'periodCount', numel(desiredOutside), 'desiredOutsideProbability', mean(desiredOutside, 'all'), ...
            'interfererInsideProbability', mean(interfererInside, 'all'), ...
            'jointPriorValidProbability', mean(~desiredOutside & ~interfererInside, 'all'));
        priorRecords = [priorRecords; row]; %#ok<AGROW>
    end
end
% 每种L都重放最差周期；保持同一快拍、先验及原型设置，不额外加载。
worstRecords = struct([]);
target = find(r.methodNames == "refined_narrow_sector");
nominalOracle = find(r.methodNames == "true_incm_nominal");
for snapshotIndex = 1:numel(s.validationSnapshotCounts)
    values = r.sinrLinear(:, :, :, snapshotIndex, target);
    [~, index] = min(values(:));
    [trial, k, point] = ind2sub(size(values), index);
    local = cfg;
    local.attitude.stdDeg(:) = s.validationJitterDeg(point);
    scene = uav.scenario.generate(local, uav.util.makeStream(cfg.random.baseSeed, s.validationSeedOffset+trial));
    obs = uav.scenario.generatePeriod(local, scene, k, uav.util.makeStream( ...
        cfg.random.baseSeed, s.validationSeedOffset+s.validationTrials+(trial-1)*s.validationPeriods+k));
    L = s.validationSnapshotCounts(snapshotIndex);
    obs.X = obs.X(:, 1:L); obs.pilot = obs.pilot(:, 1:L); obs.Rx = obs.X*obs.X'/L;
    candidate = r.candidateCfg;
    candidate.timeline.L = L;
    candidate.spectrum.desiredExclusionHalfWidthDeg = s.narrowSectorHalfWidthDeg;
    candidate.literature.guLeshem2012.upaDiagnosis.refineCandidates = true;
    [Rhat, d] = uav.incm.upaCovfit(uav.experiment.algorithmObservation(obs), candidate);
    w = uav.beamforming.solve("mvdr", Rhat, cfg.derived.nominalSteeringReference, cfg);
    m = uav.metrics.evaluate(w, obs.truth, cfg);
    reference = r.sinrLinear(trial,k,point,snapshotIndex,target);
    replayError = abs(m.sinrLinear-reference)/max(reference, eps);
    assert(replayError <= sqrt(eps), '最差周期无法从配置及子流复现。');
    outside = uav.spectrum.admissibleDirections(obs.truth.directionCosines.', candidate);
    % 仅作Oracle身份干预：不重估方向与功率，只放回被角域标签误删的干扰项。
    % 这能把“标签删除错误”与拟合、求解或名义参考失配进一步分开。
    restoredSinrDb = NaN; restoredPower = NaN; restoredSimilarity = NaN;
    interferersInside = find(~outside(2:end));
    desiredCount = size(d.desiredCandidateNodes,1);
    if ~isempty(interferersInside) && desiredCount>0
        desiredManifold = d.candidateManifold(:,1:desiredCount);
        similarity = abs(desiredManifold'*obs.truth.trueASVs(:,interferersInside+1)).^2/cfg.derived.M^2;
        [restoredSimilarity,linearIndex] = max(similarity(:));
        [candidateIndex,~] = ind2sub(size(similarity),linearIndex);
        column = desiredManifold(:,candidateIndex);
        restoredPower = d.fittedPowers(candidateIndex);
        restoredCovariance = Rhat+restoredPower*(column*column');
        restoredWeight = uav.beamforming.solve("mvdr",restoredCovariance,cfg.derived.nominalSteeringReference,cfg);
        restoredMetrics = uav.metrics.evaluate(restoredWeight,obs.truth,cfg);
        restoredSinrDb = restoredMetrics.sinrDb;
    end
    row = struct('L', L, 'trial', trial, 'period', k, 'jitterDeg', s.validationJitterDeg(point), ...
        'sinrDb', m.sinrDb, 'trueIncmNominalSinrDb', ...
        10*log10(r.sinrLinear(trial,k,point,snapshotIndex,nominalOracle)), ...
        'desiredGain', m.desiredGain, 'interferenceLeakage', m.interferenceLeakage, ...
        'relativeReplayError', replayError, 'desiredOutsidePrior', outside(1), ...
        'interfererInsidePrior', any(~outside(2:end)), 'fitResidual', d.fitResidual, ...
        'oracleRelabelSinrDb',restoredSinrDb,'restoredCandidatePower',restoredPower, ...
        'restoredCandidateSimilarity',restoredSimilarity, ...
        'usedTruth', d.usedTruth);
    worstRecords = [worstRecords; row]; %#ok<AGROW>
end
% 按预先定义的三个性能指标汇总，保留每个抖动点和每个有效鲁棒基线。
rankingRecords = struct([]);
baselineNames = ["esp", "wco", "dl_mvdr", "upa_covfit_incm"];
for snapshotIndex = 1:numel(s.validationSnapshotCounts)
    L = s.validationSnapshotCounts(snapshotIndex);
    targetRows = r.summary(r.summary.L==L & r.summary.method=="refined_narrow_sector",:);
    medianGaps = []; tailGaps = []; outageGaps = [];
    for point = 1:height(targetRows)
        baselineRows = r.summary(r.summary.L==L & ...
            r.summary.jitterDeg==targetRows.jitterDeg(point) & ismember(r.summary.method,baselineNames),:);
        medianGaps = [medianGaps;targetRows.medianSinrDb(point)-baselineRows.medianSinrDb]; %#ok<AGROW>
        tailGaps = [tailGaps;targetRows.periodLowSinrDb(point)-baselineRows.periodLowSinrDb]; %#ok<AGROW>
        outageGaps = [outageGaps;baselineRows.outageProbability-targetRows.outageProbability(point)]; %#ok<AGROW>
    end
    runtime = r.runtimeSec(:,:,:,snapshotIndex,target);
    row = struct('L',L,'minimumMedianSinrDb',min(targetRows.medianSinrDb), ...
        'maximumMedianSinrDb',max(targetRows.medianSinrDb), ...
        'minimumPeriodP05Db',min(targetRows.periodLowSinrDb), ...
        'maximumOutageProbability',max(targetRows.outageProbability), ...
        'minimumMedianGainDb',min(medianGaps),'minimumP05GainDb',min(tailGaps), ...
        'minimumOutageImprovement',min(outageGaps), ...
        'meanRuntimeMs',mean(runtime,'all')*1e3, ...
        'p95RuntimeMs',prctile(runtime(:),100*cfg.experiments.aggregation.runtimeQuantile)*1e3);
    rankingRecords = [rankingRecords;row]; %#ok<AGROW>
end
audit = struct('sourceFile', resultFile, 'dimensionCheckPassed', true, ...
    'sourceSnapshotConsistencyPassed',true, ...
    'finiteMetricsPassed', true, 'oracleBoundPassed', true, ...
    'priorStatistics', struct2table(priorRecords), 'worstCases', struct2table(worstRecords), ...
    'pairedComparisons', r.pairedComparisons, 'rankingSummary',struct2table(rankingRecords), ...
    'generatedAt', datetime('now'));
directory = fileparts(resultFile);
save(fullfile(directory, "completion_audit_"+cfg.profile+".mat"), 'audit');
writetable(audit.priorStatistics, fullfile(directory, "prior_audit_"+cfg.profile+".csv"));
writetable(audit.worstCases, fullfile(directory, "worst_cases_"+cfg.profile+".csv"));
writetable(audit.rankingSummary, fullfile(directory, "ranking_summary_"+cfg.profile+".csv"));
disp(audit.priorStatistics); disp(audit.worstCases);
disp(audit.rankingSummary);
fprintf('VALIDATION_AUDIT finite=PASS oracle=PASS reproduction=PASS\n');
end
