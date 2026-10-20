function analysis = analyzeLinYueSnapshotStudy(matPath)
%ANALYZELINYUESNAPSHOTSTUDY 配对效应区间、候选覆盖和失败率的只读结果分析。
%   不重选候选、不重新拟合、不删失败试次；区间是逐对探索性区间，未作族校正。
data=load(matPath,'result');r=data.result;cfg=r.cfg;
T=size(r.raw.sinrLinear,1);K=size(r.raw.sinrLinear,2);
contrasts=[5,6;5,10;5,19;6,20;21,22;5,22;5,16;22,23;22,24];
rows=cell(numel(r.scanPoints)*numel(r.jitterStdDeg)*size(contrasts,1),1);index=0;
for ji=1:numel(r.jitterStdDeg)
    for ip=1:numel(r.scanPoints)
        stream=uav.util.makeStream(cfg.random.baseSeed, ...
            uav.experiment.substream(cfg,7,r.originalPointIndices(ip),0,0,'bootstrap'));
        resample=randi(stream,T,T,cfg.experiments.bootstrapReplicates);
        for ci=1:size(contrasts,1)
            base=contrasts(ci,1);comparison=contrasts(ci,2);
            a=reshape(r.raw.sinrLinear(:,:,ip,ji,comparison),T,K);
            b=reshape(r.raw.sinrLinear(:,:,ip,ji,base),T,K);
            valid=isfinite(a)&isfinite(b);effect=NaN;limits=[NaN,NaN];
            if any(valid,'all')
                effect=10*log10(median(a(valid))/median(b(valid)));
                boot=nan(1,size(resample,2));
                for q=1:numel(boot)
                    aa=a(resample(:,q),:);bb=b(resample(:,q),:);ok=isfinite(aa)&isfinite(bb);
                    if any(ok,'all'),boot(q)=10*log10(median(aa(ok))/median(bb(ok)));end
                end
                limits=quantile(boot(isfinite(boot)),cfg.experiments.aggregation.confidenceQuantiles);
            end
            index=index+1;
            rows{index}=table(r.jitterStdDeg(ji),r.scanPoints(ip),r.methodNames(base),r.methodNames(comparison), ...
                effect,limits(1),limits(2),sum(valid,'all'),T*K-sum(valid,'all'), ...
                'VariableNames',{'jitterStdDeg','snapshotCount','baseline','comparison', ...
                'pairedMedianGainDb','ciLowDb','ciHighDb','validPairs','excludedFailurePairs'});
        end
    end
end
paired=vertcat(rows{:});
diagnosticRows=cell(numel(r.scanPoints)*numel(r.jitterStdDeg),1);index=0;
for ji=1:numel(r.jitterStdDeg)
    for ip=1:numel(r.scanPoints)
        index=index+1;record=struct('jitterStdDeg',r.jitterStdDeg(ji),'snapshotCount',r.scanPoints(ip));
        for field=reshape(string(fieldnames(r.diagnostics)),1,[])
            values=reshape(r.diagnostics.(field)(:,:,ip,ji),[],1);values=values(isfinite(values));
            value=NaN;if ~isempty(values),value=median(values);end
            record.(field)=value;
        end
        diagnosticRows{index}=struct2table(record);
    end
end
diagnostics=vertcat(diagnosticRows{:});
prefix="lin_yue_"+cfg.profile;
writetable(paired,fullfile(cfg.output.tableDir,prefix+"_paired_effects.csv"));
writetable(diagnostics,fullfile(cfg.output.tableDir,prefix+"_candidate_diagnostics.csv"));
failures=nnz(r.raw.status==-1);notApplicable=nnz(r.raw.status==0);
replay=r.diagnostics.replayAbsoluteError;replay=replay(isfinite(replay));
analysis=struct('paired',paired,'candidateDiagnostics',diagnostics, ...
    'failureCount',failures,'notApplicableCount',notApplicable,'maximumReplayError',max(replay), ...
    'inferenceScope',"逐对探索性配对试次簇bootstrap，非全方法同时区间", ...
    'sourceMat',matPath,'sourceCode',fileread(mfilename('fullpath')+".m"));
save(fullfile(cfg.output.rawDir,prefix+"_analysis.mat"),'analysis');
fprintf('分析完成：最大复放误差=%g，失败=%d，不适用/未计算=%d\n',analysis.maximumReplayError,failures,notApplicable);
end
