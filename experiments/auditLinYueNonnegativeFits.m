function audit = auditLinYueNonnegativeFits(matPath)
%AUDITLINYUENONNEGATIVEFITS 确定性复放全部研究样本，独立核查NNLS的KKT残差。
%   不改求解器、容差或候选，不因警告删除样本，也不重跑经典基线。
%   H=|A^H A|.^2+lambda*I，b=diag(A^H*(Rx-noise*I)*A)。
%   非负二次规划要求p>=0、g=Hp-b>=0、p.*g=0。
%   使用相对尺度衡量违反程度；即使lsqnonneg发出警告，也检查其返回点。
loaded=load(matPath,'result');r=loaded.result;cfg=r.cfg;
request=config(cfg.profile);
cfg.experiments.linYueStudy.auditKktTolerance=request.experiments.linYueStudy.auditKktTolerance;
names=["current_covfit","same_nnls0","low_load_nnls0","refined_nnls0"];
T=size(r.raw.status,1);K=size(r.raw.status,2);
records=cell(T*K*numel(r.scanPoints)*numel(r.jitterStdDeg)*numel(names),1);row=0;
for ji=1:numel(r.jitterStdDeg)
    for ip=1:numel(r.scanPoints)
        c=cfg;c.attitude.stdDeg(:)=r.jitterStdDeg(ji);c.timeline.L=r.scanPoints(ip);
        c=uav.config.validate(c);seed=r.seedScheme{ip,ji};
        for tr=1:T
            scene=uav.scenario.generate(c,uav.util.makeStream(c.random.baseSeed,seed.scenarioSubstreams(tr)));
            for k=1:K
                obs=uav.scenario.generatePeriod(c,scene,k, ...
                    uav.util.makeStream(c.random.baseSeed,seed.snapshotSubstreams(tr,k)));
                safe=uav.experiment.algorithmObservation(obs);
                [~,original]=uav.incm.upaCovfit(safe,c);
                low=c;low.incm.staticCovfit.loadingRatio=c.experiments.linYueStudy.lowLoadingRatio;
                [~,lowD]=uav.incm.upaCovfit(safe,low);
                [refinedA,~,~]=uav.spectrum.refineInterfererCandidates(safe,lowD,low);
                for q=1:numel(names)
                    A=original.candidateManifold;lambda=0;
                    if q==1
                        lambda=c.incm.staticCovfit.fitRegularization;p=original.fittedPowers;
                    else
                        if q==3,A=lowD.candidateManifold;elseif q==4,A=refinedA;end
                        p=uav.incm.fitNonnegativePowers(safe.Rx,A,original.noisePower,lambda);
                    end
                    H=abs(A'*A).^2+lambda*eye(size(A,2));
                    target=safe.Rx-original.noisePower*eye(c.derived.M);
                    b=real(sum(conj(A).*(target*A),1)).';g=H*p-b;
                    scale=max([norm(H,inf)*norm(p,inf),norm(b,inf),eps]);
                    dual=max([0;-g])/scale;
                    complement=norm(p.*g,inf)/(max(norm(p,inf),eps)*scale);
                    primal=max([0;-p])/max(norm(p,inf),eps);
                    residual=max([dual,complement,primal]);
                    row=row+1;
                    records{row}=table(r.jitterStdDeg(ji),r.scanPoints(ip),tr,k,names(q), ...
                        residual,dual,complement,primal,cond(H), ...
                        'VariableNames',{'jitterStdDeg','snapshotCount','trial','period','method', ...
                        'relativeKktResidual','dualViolation','complementarity','primalViolation','gramCondition'});
                end
            end
        end
        fprintf('NNLS审计：抖动%g°，L=%d完成\n',r.jitterStdDeg(ji),r.scanPoints(ip));
    end
end
detail=vertcat(records{:});summaries=cell(numel(names),1);
for q=1:numel(names)
    value=detail.relativeKktResidual(detail.method==names(q));
    summaries{q}=table(names(q),numel(value),median(value),max(value), ...
        sum(value>cfg.experiments.linYueStudy.auditKktTolerance), ...
        'VariableNames',{'method','count','medianKktResidual','maxKktResidual','exceedsTolerance'});
end
summary=vertcat(summaries{:});prefix="lin_yue_"+cfg.profile+"_nnls_audit";
writetable(detail,fullfile(cfg.output.tableDir,prefix+"_detail.csv"));
writetable(summary,fullfile(cfg.output.tableDir,prefix+"_summary.csv"));
audit=struct('detail',detail,'summary',summary,'sourceMat',matPath, ...
    'tolerance',cfg.experiments.linYueStudy.auditKktTolerance, ...
    'sourceCode',fileread(mfilename('fullpath')+".m"));
save(fullfile(cfg.output.rawDir,prefix+".mat"),'audit');disp(summary);
end
