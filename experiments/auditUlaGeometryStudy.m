function audit = auditUlaGeometryStudy(profile)
%AUDITULAGEOMETRYSTUDY 检查一维数据完整性，并导出同世界源的阵型诊断。
arguments
    profile (1,1) string = "preliminary"
end
cfg=config(profile); g=cfg.experiments.recUlaGeometry;
root=fileparts(fileparts(mfilename('fullpath')));
directory=fullfile(root,g.outputDirectory,profile);
loaded=load(fullfile(directory,'ula_geometry.mat'),'results'); r=loaded.results; t=r.records;
expected=numel(r.geometries)*numel(g.halfWidthsDeg)*(numel(g.rawGridSizes)+2+ ...
    g.trialCount*numel(g.snapshotCounts)*3);
assert(height(t)==expected,'ULA记录数不完整。');
invalid=t.method=="raw_rec" & t.L<g.sensorCount;
assert(all(isnan(t.sinrLinear(invalid))) && all(t.status(invalid)=="not-applicable-L<M"));
assert(all(isfinite(t.sinrLinear(~invalid))) && all(t.oracleGapDb(~invalid)>=-1e-7));
assert(all(t.weightLoading(t.method=="raw_rec" & ~invalid)==0),'原式触发额外加载。');
convergence=struct([]);
for index=1:numel(r.geometries)
    for half=g.halfWidthsDeg
        rows=t(t.geometryIndex==index & isinf(t.L) & t.method=="raw_rec" & t.halfWidthDeg==half,:);
        entry=struct('geometryIndex',index,'halfWidthDeg',half, ...
            'coarseSinrDb',rows.sinrDb(1),'fineSinrDb',rows.sinrDb(end), ...
            'differenceDb',rows.sinrDb(end)-rows.sinrDb(1));
        convergence=[convergence;entry]; %#ok<AGROW>
    end
end
writetable(struct2table(convergence),fullfile(directory,'ula_grid_comparison.csv'));
% 相同世界源、阵元数及功率，但孔径、可观测维数、先验含义随阵型改变。
% 只作诊断表，不能用其均值证明“维度导致性能更差”。
previous=load(fullfile(root,cfg.experiments.recGeometry.outputDirectory,profile,'geometry_study.mat'),'results');
u=previous.results.records; comparisons=struct([]);
matched=find(string({r.geometries.family})=="matched_world_projection");
matchedRows=t.family=="matched_world_projection";
assert(all(t.snapshotSubstream(matchedRows)==cfg.experiments.recGeometry.seedOffset+t.trial(matchedRows)+1), ...
    '同世界源对照没有使用UPA配对子流。');
for index=matched
    worldIndex=r.geometries(index).parameter;
    for L=[Inf,g.snapshotCounts]
        for half=g.halfWidthsDeg
            for pair=["raw_rec","covfit_refined"]
                ul=t.geometryIndex==index & t.L==L & t.halfWidthDeg==half & t.method==pair;
                if pair=="raw_rec"
                    upMethod="raw_rec_"+half;
                    if isinf(L), ul=ul & t.rawGridSize==max(g.rawGridSizes); end
                else
                    upMethod="refined_covfit_"+half;
                end
                up=u.geometryIndex==worldIndex & u.L==L & u.method==upMethod;
                up=up & u.trial<=g.trialCount;
                if pair=="raw_rec" && isinf(L)
                    up=up & u.rawGridSize==max(cfg.experiments.recGeometry.rawGridSizes);
                end
                uOracle=u.sinrLinear(u.geometryIndex==worldIndex & u.L==L & u.method=="oracle" & u.trial<=g.trialCount);
                entry=struct('worldGeometryIndex',worldIndex,'L',L,'halfWidthDeg',half,'method',pair, ...
                    'ulaSinrDb',10*log10(median(t.sinrLinear(ul))), ...
                    'upaSinrDb',10*log10(median(u.sinrLinear(up))), ...
                    'ulaOracleDb',median(t.oracleSinrDb(ul)), ...
                    'upaOracleDb',10*log10(median(uOracle)));
                entry.ulaOracleGapDb=entry.ulaOracleDb-entry.ulaSinrDb;
                entry.upaOracleGapDb=entry.upaOracleDb-entry.upaSinrDb;
                comparisons=[comparisons;entry]; %#ok<AGROW>
            end
        end
    end
end
writetable(struct2table(comparisons),fullfile(directory,'matched_world_comparison.csv'));
% 对预设近邻干扰族作同迹拆分干预：只把已拟合的总干扰功率重新放到
% 两个真实干扰指纹上，噪声和最终期望参考不变。此项读取真值，仅诊断。
splitRecords=struct([]);
pairIndices=find(string({r.geometries.family})=="interferer_pair_separation");
for index=pairIndices
    A=uav.array.steeringULA(deg2rad(r.geometries(index).anglesDeg),g.sensorCount,g.spacingOverLambda);
    Ri=r.cfg.derived.noisePower*eye(g.sensorCount)+A(:,2:end)*diag(r.cfg.derived.interfererPowers)*A(:,2:end)';
    for halfIndex=1:numel(g.halfWidthsDeg)
        d=r.candidateDiagnostics{index,halfIndex,2};
        mass=sum(d.powers(d.desiredCount+1:end));
        assigned=mass*r.cfg.derived.interfererPowers(:)/sum(r.cfg.derived.interfererPowers);
        Rh=A(:,2:end)*diag(assigned)*A(:,2:end)'+d.noisePower*eye(g.sensorCount);
        w=uav.beamforming.solve('mvdr',Rh,A(:,1),r.cfg);
        select=t.geometryIndex==index & isinf(t.L) & t.method=="covfit_refined" & ...
            t.halfWidthDeg==g.halfWidthsDeg(halfIndex);
        entry=struct('geometryIndex',index,'separationDeg',r.geometries(index).parameter, ...
            'halfWidthDeg',g.halfWidthsDeg(halfIndex),'fittedInterfererPowerSum',mass, ...
            'originalSinrDb',t.sinrDb(select), ...
            'truthSplitSinrDb',10*log10(r.cfg.derived.desiredPower*abs(w'*A(:,1))^2/real(w'*Ri*w)), ...
            'interfererInsideSector',t.interfererInsideSector(select), ...
            'usedTruthForDiagnosisOnly',true);
        splitRecords=[splitRecords;entry]; %#ok<AGROW>
    end
end
writetable(struct2table(splitRecords),fullfile(directory,'ula_truth_split_diagnosis.csv'));
audit=struct('passed',true,'recordCount',height(t),'geometryCount',numel(r.geometries), ...
    'notApplicableCount',sum(invalid),'maxGridDifferenceDb',max(abs([convergence.differenceDb])), ...
    'generatedAt',datetime('now'));
save(fullfile(directory,'ula_audit.mat'),'audit'); disp(audit);
end
