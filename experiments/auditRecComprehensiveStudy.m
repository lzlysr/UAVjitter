function audit=auditRecComprehensiveStudy(directory)
%AUDITRECCOMPREHENSIVESTUDY 独立复核落盘完整性、约束、状态与跨条件配对。
load(fullfile(directory,'manifest.mat'),'manifest');g=manifest.cfg.experiments.recComprehensive;
results=cell(1,numel(manifest.cases)); audit=struct('caseCount',numel(results), ...
    'observations',0,'successfulMethodEvaluations',0,'notApplicable',0,'numericalFailures',0, ...
    'kktViolations',0,'pairedConditionChecks',0,'replayChecks',0);
for id=1:numel(results)
    path=fullfile(directory,sprintf('case_%06d.mat',id));assert(isfile(path),'有条件块未完成。');
    data=load(path,'result');r=data.result;results{id}=r;c=r.condition;
    assert(strcmp(r.sourceSignature,manifest.sourceSignature));
    assert(isequal(c,manifest.cases(id)));
    ok=r.status==1; audit.observations=audit.observations+c.trials*numel(c.snapshotCounts);
    audit.successfulMethodEvaluations=audit.successfulMethodEvaluations+sum(ok,'all');
    audit.notApplicable=audit.notApplicable+sum(r.status==0,'all');
    audit.numericalFailures=audit.numericalFailures+sum(r.status==-1,'all');
    audit.kktViolations=audit.kktViolations+sum(r.fitKkt(ok)>g.kktTolerance);
    assert(all(isfinite(r.sinrLinear(ok)))&&all(r.sinrLinear(ok)>0));
    assert(all(r.oracleLossDb(ok)>=-1e-6));
    if c.jitterStdDeg==0,assert(all(abs(r.desiredGain(ok)-1)<1e-7));end
    expected=c.trials*sum(c.snapshotCounts<g.sensorCount)*sum(r.methodNames=="raw_rec");
    assert(sum(r.status==0,'all')==expected,'不适用状态与预定的L<M条件不符。');
end
for ai=1:numel(g.arrayTypes)
    ids=find(arrayfun(@(c)c.block=="snapshots"&&c.arrayType==g.arrayTypes(ai),manifest.cases));
    for id=ids
        base=results{id}; c=base.condition;
        other=find(arrayfun(@(q)q.geometryIndex==c.geometryIndex&&q.arrayType==c.arrayType&& ...
            (q.block=="power"||q.block=="sector"||q.block=="jitter")&&q.snrDb==g.baseSnrDb&& ...
            isequal(q.inrDb,g.baseInrDb)&&q.halfWidthDeg==g.baseHalfWidthDeg&&q.jitterStdDeg==0,manifest.cases));
        for j=other
            r=results{j};[~,positions]=ismember(r.condition.snapshotCounts,c.snapshotCounts);
            assert(isequaln(r.sinrLinear,base.sinrLinear(1:r.condition.trials,positions,:)), ...
                '同条件跨区块未能逐元素配对。');
            assert(isequal(r.snapshotSubstreams,base.snapshotSubstreams(1:r.condition.trials)));
            audit.pairedConditionChecks=audit.pairedConditionChecks+1;
        end
    end
    % 预先固定复跑每阵型第一个主快拍条件的前两个试次，逐元素检查数据。
    if isempty(ids)
        ids=find(arrayfun(@(c)c.arrayType==g.arrayTypes(ai),manifest.cases),1);
    end
    if isempty(ids),continue;end
    id=ids(1);c=manifest.cases(id);c.trials=min(2,c.trials);
    repeated=uav.experiment.runRecComprehensiveCase(manifest.cfg,manifest.geometries(c.geometryIndex),c);
    assert(isequaln(repeated.sinrLinear,results{id}.sinrLinear(1:c.trials,:,:)));
    audit.replayChecks=audit.replayChecks+1;
end
audit.sourceSignature=manifest.sourceSignature;
audit.postprocessingSource=fileread(mfilename('fullpath')+".m");
audit.complete=true;save(fullfile(directory,'audit.mat'),'audit');disp(audit);
end
