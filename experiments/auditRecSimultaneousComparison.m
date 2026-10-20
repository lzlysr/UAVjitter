function audit = auditRecSimultaneousComparison(resultFile)
%AUDITRECSIMULTANEOUSCOMPARISON 用同一试次重采样核查整组中位数比较。
%   比较族固定为细化7°候选对ESP/WCO/DL/历史CovFit的全部L×抖动点。
%   以bootstrap标准差标准化误差，再取跨比较最大绝对误差的95%分位，
%   构造近似同时区间。它是有限样本bootstrap审计，不是精确覆盖率定理。
arguments
    resultFile (1,1) string
end
saved = load(resultFile, 'results');
r = saved.results;
cfg = r.cfg;
s = cfg.experiments.recRootCause;
target = find(r.methodNames=="refined_narrow_sector");
references = ["esp", "wco", "dl_mvdr", "upa_covfit_incm"];
T = size(r.sinrLinear,1);
B = cfg.experiments.bootstrapReplicates;
stream = uav.util.makeStream(cfg.random.baseSeed,s.validationSeedOffset*4);
indices = randi(stream,T,T,B);
pointCount = numel(s.validationSnapshotCounts)*numel(s.validationJitterDeg)*numel(references);
bootDifferences = zeros(pointCount,B);
records = struct([]);
point = 0;
for snapshotIndex = 1:numel(s.validationSnapshotCounts)
    for jitterIndex = 1:numel(s.validationJitterDeg)
        targetMean = mean(r.sinrLinear(:,:,jitterIndex,snapshotIndex,target),2);
        targetBootstrap = 10*log10(median(targetMean(indices),1));
        for reference = references
            point = point+1;
            referenceMean = mean(r.sinrLinear(:,:,jitterIndex,snapshotIndex,r.methodNames==reference),2);
            bootDifferences(point,:) = targetBootstrap-10*log10(median(referenceMean(indices),1));
            row = struct('L',s.validationSnapshotCounts(snapshotIndex), ...
                'jitterDeg',s.validationJitterDeg(jitterIndex),'reference',reference, ...
                'differenceDb',10*log10(median(targetMean)/median(referenceMean)));
            records = [records;row]; %#ok<AGROW>
        end
    end
end
output = struct2table(records);
bootstrapSd = std(bootDifferences,0,2);
errors = bootDifferences-output.differenceDb;
standardized = zeros(size(errors));
nonzero = bootstrapSd>eps(max(abs(output.differenceDb),1));
standardized(nonzero,:) = errors(nonzero,:)./bootstrapSd(nonzero);
confidence = 1-2*cfg.experiments.aggregation.confidenceQuantiles(1);
critical = quantile(max(abs(standardized),[],1),confidence);
output.bootstrapSdDb = bootstrapSd;
output.simultaneousLowerDb = output.differenceDb-critical*bootstrapSd;
output.simultaneousUpperDb = output.differenceDb+critical*bootstrapSd;
audit = struct('target',r.methodNames(target),'references',references, ...
    'comparisonCount',pointCount,'independentTrialCount',T,'bootstrapReplicates',B, ...
    'confidenceLevel',confidence,'criticalMaximumStandardizedError',critical, ...
    'records',output,'bootstrapDifferences',bootDifferences, ...
    'finiteSampleApproximation',true,'generatedAt',datetime('now'));
directory = fileparts(resultFile);
save(fullfile(directory,"simultaneous_comparisons_"+cfg.profile+".mat"),'audit');
writetable(output,fullfile(directory,"simultaneous_comparisons_"+cfg.profile+".csv"));
fprintf('SIMULTANEOUS_AUDIT comparisons=%d positiveLower=%d minLower=%.6g\n', ...
    pointCount,sum(output.simultaneousLowerDb>0),min(output.simultaneousLowerDb));
disp(output(output.simultaneousLowerDb<=0,:));
end
