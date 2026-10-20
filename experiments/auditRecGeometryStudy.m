function audit = auditRecGeometryStudy(profile)
%AUDITRECGEOMETRYSTUDY 检查几何扫描的完整性、重复基准及原式适用性。
arguments
    profile (1,1) string = "smoke"
end
cfg = config(profile);
g = cfg.experiments.recGeometry;
root = fileparts(fileparts(mfilename('fullpath')));
directory = fullfile(root,g.outputDirectory,profile);
loaded = load(fullfile(directory,'geometry_study.mat'),'results');
r = loaded.results;
t = r.records;
geometryCount = numel(r.geometries);
% 总体每点7种方法，其中两种Rec各多一个加密网格；有限样本每L均7条记录。
expectedCount = geometryCount*(numel(r.methodNames)+numel(g.rawGridSizes)-1+ ...
    numel(g.rawGridSizes)-1+g.trialCount*numel(g.snapshotCounts)*numel(r.methodNames));
assert(height(t)==expectedCount,'几何扫描记录数不完整。');
raw = startsWith(t.method,"raw_rec");
inapplicable = raw & t.L<r.cfg.derived.M;
assert(all(t.status(inapplicable)=="not-applicable-L<M"),'欠采样原式状态错误。');
assert(all(isnan(t.sinrLinear(inapplicable))),'不得暗中补齐欠采样原式。');
assert(all(isfinite(t.sinrLinear(~inapplicable))),'可用方法出现非有限值。');
assert(all(t.weightLoading(raw & ~inapplicable)==0),'原式权值触发了额外数值加载，需单独披露。');
% 四个重复默认点的几何和数据应一致；不检查运行时间等非确定性量。
defaultIndices = find(arrayfun(@(p) isequal(p.desired,r.cfg.sources.desiredAzElDeg) && ...
    isequal(p.interferers,r.cfg.sources.interfererAzElDeg),r.geometries));
reference = t(t.geometryIndex==defaultIndices(1),:);
maxDuplicateDifference = 0;
for index = defaultIndices(2:end)
    duplicate = t(t.geometryIndex==index,:);
    assert(isequal(reference.method,duplicate.method) && isequal(reference.trial,duplicate.trial));
    delta = abs(reference.sinrLinear-duplicate.sinrLinear);
    maxDuplicateDifference = max(maxDuplicateDifference,max(delta,[],'omitnan'));
end
assert(maxDuplicateDifference==0,'重复基准未实现逐元素复现。');
% 另存总体候选诊断：真值只用于事后匹配，不进入候选生成或功率拟合。
candidateAudit = cell(geometryCount,3);
candidateMethods = ["historical_covfit","refined_covfit_25","refined_covfit_7"];
study = r.cfg.experiments.recRootCause;
for geometryIndex = 1:geometryCount
    local = r.cfg;
    local.sources.desiredAzElDeg = r.geometries(geometryIndex).desired;
    local.sources.interfererAzElDeg = r.geometries(geometryIndex).interferers;
    local = uav.config.validate(local);
    scenario = uav.scenario.generate(local,uav.util.makeStream(local.random.baseSeed,g.seedOffset));
    obs = uav.scenario.generatePeriod(local,scenario,1,uav.util.makeStream(local.random.baseSeed,g.seedOffset+1));
    A = obs.truth.trueASVs;
    obs.Rx = obs.truth.trueIncm + local.derived.desiredPower*A(:,1)*A(:,1)';
    safe = uav.experiment.algorithmObservation(obs);
    for methodIndex = 1:numel(candidateMethods)
        point = local;
        if methodIndex>1
            point.spectrum.gridSize = study.candidateCovfit.gridSize;
            point.literature.guLeshem2012.upaDiagnosis.loadingRatio = study.candidateCovfit.loadingRatio;
            point.literature.guLeshem2012.upaDiagnosis.peakNmsRadiusBins = study.candidateCovfit.peakNmsRadiusBins;
            point.literature.guLeshem2012.upaDiagnosis.fitRegularization = study.candidateCovfit.fitRegularization;
            point.literature.guLeshem2012.upaDiagnosis.refineCandidates = true;
            if methodIndex==3
                point.spectrum.desiredExclusionHalfWidthDeg = study.narrowSectorHalfWidthDeg;
            end
            point = uav.config.validate(point);
        end
        [Rh,d] = uav.incm.upaCovfit(safe,point);
        fingerprint = abs(A'*d.candidateManifold).^2/point.derived.M^2;
        [similarity,sourceIndex] = max(fingerprint,[],1);
        outside = uav.spectrum.admissibleDirections(point.derived.nominalDirectionCosines.',point);
        candidateAudit{geometryIndex,methodIndex} = struct('method',candidateMethods(methodIndex), ...
            'diagnostics',d,'reconstructedIncm',Rh,'truthFingerprintSimilarity',fingerprint, ...
            'nearestTrueSourceIndex',sourceIndex,'nearestSimilarity',similarity, ...
            'trueInterfererInsideDesiredSector',~outside(2:end));
    end
end
save(fullfile(directory,'population_candidate_audit.mat'),'candidateAudit','candidateMethods','-v7.3');
audit = struct('profile',profile,'recordCount',height(t),'geometryCount',geometryCount, ...
    'inapplicableCount',sum(inapplicable),'defaultGeometryIndices',defaultIndices, ...
    'maxDuplicateDifference',maxDuplicateDifference,'passed',true,'generatedAt',datetime('now'));
save(fullfile(directory,'geometry_audit.mat'),'audit');
disp(audit);
end
