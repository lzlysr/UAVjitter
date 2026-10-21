function comparison = compareCovfitRefinements(matPath)
%COMPARECOVFITREFINEMENTS 对同批观测比较两种Le移植与本次局部细化。
%   同时给出配对中位数比、试次簇bootstrap、逐周期差异和失败数。
%   0.1dB只是事先明确的数值接近描述阈值，不冒称普适统计等效性。
loaded=load(matPath,'result');r=loaded.result;cfg=r.cfg;
% 兼容在新增描述阈值前导出的smoke；阈值只影响接近比例，不改任何仿真值。
if ~isfield(cfg.experiments.linYueStudy,'comparisonToleranceDb')
    request=config(cfg.profile);
    cfg.experiments.linYueStudy.comparisonToleranceDb=request.experiments.linYueStudy.comparisonToleranceDb;
end
reference=find(r.methodNames=="refined_nnls0",1);
names=["le_refined_covfit","le_shared_candidates"];
T=size(r.raw.sinrLinear,1);K=size(r.raw.sinrLinear,2);rows={};
parameterStudy=isfield(r,'axisNames');
if parameterStudy,blocks=numel(r.axisNames);else,blocks=numel(r.jitterStdDeg);end
for block=1:blocks
    if parameterStudy,pointCount=numel(r.scanPoints{block});else,pointCount=numel(r.scanPoints);end
    for ip=1:pointCount
        if parameterStudy
            parameter=r.axisNames(block);value=r.scanPoints{block}(ip);jitter=cfg.attitude.stdDeg(1);
        else
            parameter="snapshotCount";value=r.scanPoints(ip);jitter=r.jitterStdDeg(block);
        end
        stream=uav.util.makeStream(cfg.random.baseSeed,uav.experiment.substream(cfg,7,ip,0,0,'bootstrap'));
        draws=randi(stream,T,T,cfg.experiments.bootstrapReplicates);
        a=reshape(r.raw.sinrLinear(:,:,ip,block,reference),T,K);
        for method=names
            index=find(r.methodNames==method,1);b=reshape(r.raw.sinrLinear(:,:,ip,block,index),T,K);
            valid=isfinite(a)&isfinite(b);differenceDb=10*log10(b(valid)./a(valid));
            effect=NaN;ci=[NaN,NaN];d=[NaN,NaN,NaN,NaN];within=NaN;
            if any(valid,'all')
                effect=10*log10(median(b(valid))/median(a(valid)));
                boot=nan(1,size(draws,2));
                for q=1:numel(boot)
                    aa=a(draws(:,q),:);bb=b(draws(:,q),:);ok=isfinite(aa)&isfinite(bb);
                    if any(ok,'all'),boot(q)=10*log10(median(bb(ok))/median(aa(ok)));end
                end
                ci=quantile(boot(isfinite(boot)),cfg.experiments.aggregation.confidenceQuantiles);
                d=[median(differenceDb),quantile(differenceDb,[.05,.95]),max(abs(differenceDb))];
                within=mean(abs(differenceDb)<=cfg.experiments.linYueStudy.comparisonToleranceDb);
            end
            sa=r.raw.status(:,:,ip,block,reference);sb=r.raw.status(:,:,ip,block,index);
            rows{end+1}=table(parameter,value,jitter,method,effect,ci(1),ci(2), ...
                d(1),d(2),d(3),d(4),within,nnz(valid),nnz(sa==-1),nnz(sb==-1), ...
                'VariableNames',{'parameter','value','jitterStdDeg','comparison', ...
                'medianGainOverLocalDb','ciLowDb','ciHighDb','medianPeriodDifferenceDb', ...
                'periodDifferenceP05Db','periodDifferenceP95Db','maxAbsolutePeriodDifferenceDb', ...
                'fractionWithinTolerance','validPairs','localFailures','comparisonFailures'}); %#ok<AGROW>
        end
    end
end
comparison=vertcat(rows{:});prefix=r.experimentName+"_"+cfg.profile+"_refinement_paired";
writetable(comparison,fullfile(cfg.output.tableDir,prefix+".csv"));
sourceCode=fileread(mfilename('fullpath')+".m"); %#ok<NASGU>
save(fullfile(cfg.output.rawDir,prefix+".mat"),'comparison','sourceCode');
disp(comparison(:,{'parameter','value','jitterStdDeg','comparison','medianGainOverLocalDb','validPairs','comparisonFailures'}));
end
