function analysis=analyzeRecComprehensiveStudy(directory)
%ANALYZERECCOMPREHENSIVESTUDY 保留所有几何，分层配对推断，不按好坏删点。
%   总体输入无抽样区间；逐条件中位数区间为次序统计量区间。
load(fullfile(directory,'manifest.mat'),'manifest');
cfg=manifest.cfg; g=cfg.experiments.recComprehensive;
records=struct([]); failures=struct([]); loaded=cell(1,numel(manifest.cases)); missing=[];
for id=1:numel(manifest.cases)
    path=fullfile(directory,sprintf('case_%06d.mat',id));
    if ~isfile(path),missing(end+1)=id;continue;end %#ok<AGROW>
    data=load(path,'result'); r=data.result;
    assert(strcmp(r.sourceSignature,manifest.sourceSignature),'结果与源码签名不一致。');
    loaded{id}=r; c=r.condition;
    for li=1:numel(c.snapshotCounts)
        for mi=1:numel(r.methodNames)
            s=r.status(:,li,mi); ok=s==1; values=10*log10(r.sinrLinear(ok,li,mi));
            linearValues=r.sinrLinear(ok,li,mi);
            losses=r.oracleLossDb(ok,li,mi); ci=uav.experiment.recMedianInterval(values,g.confidenceLevel);
            if c.population && any(ok),ci=[median(values),median(values)];end
            active=sum(s~=0); outage=sum(losses>g.oracleLossThresholdDb)+sum(s==-1);
            if active>0
                rate=outage/active; z=sqrt(2)*erfcinv(1-g.confidenceLevel);
                middle=(rate+z^2/(2*active))/(1+z^2/active);
                radius=z*sqrt(rate*(1-rate)/active+z^2/(4*active^2))/(1+z^2/active);
                outageCI=[middle-radius,middle+radius];
            else,rate=NaN;outageCI=[NaN,NaN];end
            if c.population,outageCI(:)=NaN;end
            q=struct('caseId',id,'block',c.block,'arrayType',c.arrayType, ...
                'geometryId',r.geometry.id,'geometryFamily',r.geometry.family, ...
                'method',r.methodNames(mi),'L',c.snapshotCounts(li),'snrDb',c.snrDb, ...
                'inrDb',firstOrNaN(c.inrDb),'J',numel(c.inrDb),'halfWidthDeg',c.halfWidthDeg, ...
                'secondInrDb',lastOrNaN(c.inrDb), ...
                'jitterStdDeg',c.jitterStdDeg,'gridSize',c.gridSize,'population',c.population, ...
                'nominalPriorValid',r.nominalPriorValid,'maxDesiredCoherence',r.maxDesiredCoherence, ...
                'successes',sum(ok),'numericalFailures',sum(s==-1),'notApplicable',sum(s==0), ...
                'medianSinrDb',10*log10(median(linearValues)),'medianCiLowDb',ci(1),'medianCiHighDb',ci(2), ...
                'p05SinrDb',10*log10(quantileOrNaN(linearValues,.05)), ...
                'medianLossDb',-10*log10(median(linearValues./r.oracleSinrLinear(ok))), ...
                'outageRate',rate,'outageCiLow',outageCI(1),'outageCiHigh',outageCI(2), ...
                'medianIncmNmse',median(r.incmNmse(ok,li,mi)), ...
                'medianShapeNmse',median(r.shapeNmse(ok,li,mi)), ...
                'medianDesiredGain',median(r.desiredGain(ok,li,mi)), ...
                'medianInterferenceLeakage',median(r.interferenceLeakage(ok,li,mi)), ...
                'medianNoiseOutput',median(r.noiseOutput(ok,li,mi)), ...
                'medianMinEigenvalue',median(r.minimumEigenvalue(ok,li,mi)), ...
                'medianConditionNumber',median(r.conditionNumber(ok,li,mi)), ...
                'medianRuntimeSec',median(r.runtimeSec(ok,li,mi)), ...
                'weightLoadedCount',sum(r.weightLoading(ok,li,mi)>0), ...
                'kktViolations',sum(r.fitKkt(ok,li,mi)>g.kktTolerance), ...
                'maxKkt',maxOrNaN(r.fitKkt(ok,li,mi)), ...
                'desiredOutsideRate',mean(r.desiredOutside),'interfererInsideRate',mean(r.interfererInside), ...
                'oracleMedianDb',median(10*log10(r.oracleSinrLinear)), ...
                'trueIncmNominalMedianDb',median(10*log10(r.trueIncmNominalSinrLinear)));
            if isempty(records),records=q;else,records(end+1)=q;end %#ok<AGROW>
            for ti=find(s==-1).'
                e=r.errors{ti,li,mi};
                failure=struct('caseId',id,'trial',ti,'L',c.snapshotCounts(li), ...
                    'method',r.methodNames(mi),'identifier',string(e.identifier),'message',string(e.message)); %#ok<AGROW>
                if isempty(failures),failures=failure;else,failures(end+1)=failure;end %#ok<AGROW>
            end
        end
    end
end
assert(isempty(missing),'uav:experiment:IncompleteStudy', ...
    '存在未完成条件块，停止总体推断；请先续跑，不按已完成子集替代预定总体。');
summary=struct2table(records); writetable(summary,fullfile(directory,'condition_summary.csv'));
if ~isempty(failures),writetable(struct2table(failures),fullfile(directory,'numerical_failures.csv'));end
[primary,bootstrap]=primaryComparisons(loaded,cfg);
if width(primary)>0,writetable(primary,fullfile(directory,'primary_paired_comparisons.csv'));end
analysis=struct('summary',summary,'primary',primary,'bootstrap',bootstrap,'missingCases',missing, ...
    'sourceSignature',manifest.sourceSignature,'createdAt',char(datetime('now')));
analysis.postprocessingSource=fileread(mfilename('fullpath')+".m");
analysis.intervalSource=fileread(which('uav.experiment.recMedianInterval'));
save(fullfile(directory,'analysis.mat'),'analysis','-v7.3');
fprintf('已汇总%d/%d块、%d行；缺失%d块，数值失败%d项。\n', ...
    sum(~cellfun(@isempty,loaded)),numel(loaded),height(summary),numel(missing),numel(failures));
disp(primary);
end

function [tableOut,boot]=primaryComparisons(loaded,cfg)
%PRIMARYCOMPARISONS 几何为上层单位，试次为下层；方法和L共用重采样索引。
g=cfg.experiments.recComprehensive; stream=uav.util.makeStream(cfg.random.baseSeed,g.bootstrapSubstream);
rows=struct([]); boot=[]; draws=cell(1,numel(g.arrayTypes)); pairs=nchoosek(1:numel(g.methodNames),2);
% 同一世界几何用于两个阵型：整个比较族共享上层抽样和下层试次索引。
% 各阵型再按自己的先验有效标记分层，不能分别重采样而丢掉阵型间相关性。
geometryDraws=randi(stream,g.randomGeometryCount,g.bootstrapReplicates,g.randomGeometryCount);
trialDraws=randi(stream,g.trialCount,g.trialCount,g.randomGeometryCount,g.bootstrapReplicates);
for ai=1:numel(g.arrayTypes)
    subset=loaded(cellfun(@(r) ~isempty(r) && r.condition.block=="snapshots" && ...
        r.condition.arrayType==g.arrayTypes(ai) && r.nominalPriorValid,loaded));
    if isempty(subset),continue;end
    validL=find(g.snapshotCounts>=g.sensorCount); lengths=g.snapshotCounts(validL);
    cube=nan(g.trialCount,numel(validL),numel(g.methodNames),numel(subset));
    for gi=1:numel(subset)
        cube(:,:,:,gi)=subset{gi}.sinrLinear(:,validL,:)./subset{gi}.oracleSinrLinear;
    end
    % 主比较不以成功试次条件化；任意失败导致对应端点不可推断。
    perGeometry=-10*log10(median(cube,1));
    point=reshape(mean(perGeometry,4),numel(validL),numel(g.methodNames));
    methodDraws=nan(g.bootstrapReplicates,numel(validL),numel(g.methodNames));
    lookup=zeros(1,g.randomGeometryCount);
    validIds=cellfun(@(r)r.geometry.randomIndex,subset);
    lookup(validIds)=1:numel(subset);
    for bi=1:g.bootstrapReplicates
        allSelected=lookup(geometryDraws(bi,:));positions=find(allSelected>0);
        selected=allSelected(positions);
        if isempty(selected),continue;end
        sample=nan(numel(selected),numel(validL),numel(g.methodNames));
        for gj=1:numel(selected)
            trialIndices=trialDraws(:,positions(gj),bi);
            sample(gj,:,:)=-10*log10(median(cube(trialIndices,:,:,selected(gj)),1));
        end
        methodDraws(bi,:,:)=mean(sample,1);
    end
    draws{ai}=methodDraws;
    for li=1:numel(lengths)
        for pi=1:size(pairs,1)
            a=pairs(pi,1); b=pairs(pi,2);
            difference=squeeze(methodDraws(:,li,a)-methodDraws(:,li,b));
            row=struct('arrayType',g.arrayTypes(ai),'L',lengths(li), ...
                'methodA',g.methodNames(a),'methodB',g.methodNames(b),'geometryCount',numel(subset), ...
                'differenceLossDb',point(li,a)-point(li,b),'bootstrapSE',std(difference,'omitnan'), ...
                'bootstrapValidReplicates',sum(isfinite(difference)), ...
                'simultaneousLowDb',NaN,'simultaneousHighDb',NaN); %#ok<AGROW>
            if isempty(rows),rows=row;else,rows(end+1)=row;end %#ok<AGROW>
            boot(:,end+1)=difference; %#ok<AGROW>
        end
    end
end
tableOut=struct2table(rows); if isempty(rows),return;end
center=tableOut.differenceLossDb.'; scale=tableOut.bootstrapSE.';
valid=isfinite(center)&isfinite(scale)&scale>0;
if any(valid)
    standardized=abs((boot(:,valid)-center(valid))./scale(valid));
    standardized=standardized(all(isfinite(standardized),2),:);
    critical=quantile(max(standardized,[],2),g.confidenceLevel);
    tableOut.simultaneousLowDb(valid)=center(valid).'-critical*scale(valid).';
    tableOut.simultaneousHighDb(valid)=center(valid).'+critical*scale(valid).';
end
% 负差值表示A的Oracle损失更小；区间跨零不是等效证据。
end

function x=firstOrNaN(v)
if isempty(v),x=NaN;else,x=v(1);end
end
function x=lastOrNaN(v)
if numel(v)<2,x=NaN;else,x=v(2);end
end
function x=quantileOrNaN(v,p)
if isempty(v),x=NaN;else,x=quantile(v,p);end
end
function x=maxOrNaN(v)
v=v(isfinite(v));if isempty(v),x=NaN;else,x=max(v);end
end
