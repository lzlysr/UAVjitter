function result = runLinYueSnapshotStudy(profile)
%RUNLINYUESNAPSHOTSTUDY 从目标图的MAT配置复放，受控比较Lin/Yue重构核心。
%   不修正最终期望导向；真值仅用于Oracle诊断及指标评价。
%   论文完整候选前端、极化阵列及Lin旁瓣约束均不冒充已复现。
arguments
    profile (1,1) string = "smoke"
end
timer=tic;
request=config(profile); g=request.experiments.linYueStudy;
source=load(g.sourceMat,'result'); archived=source.result;
cfg=archived.cfg;
assert(cfg.array.Nx==3 && cfg.array.Ny==3,'目标MAT必须是3×3 UPA。');
assert(isequal(cfg.attitude.stdDeg,[5,5,5]),'当前协议复放三轴5°目标图。');
cfg.experiments.linYueStudy=g;
cfg.profile=profile;
if profile=="smoke"
    cfg.monteCarlo.nTrials=min(g.smokeTrials,cfg.monteCarlo.nTrials);
    cfg.experiments.periodCounts.snapshotSweep=min(g.smokePeriods,cfg.experiments.periodCounts.snapshotSweep);
elseif profile~="preliminary"
    error('uav:study:Profile','本研究只提供smoke与目标MAT规模的preliminary。');
end
root=fileparts(fileparts(mfilename('fullpath')));
cfg.output=struct('rootDir',root,'overwrite',true, ...
    'rawDir',fullfile(root,'results','raw',g.outputTag), ...
    'tableDir',fullfile(root,'results','tables',g.outputTag), ...
    'figureDir',fullfile(root,'results','figures',g.outputTag), ...
    'logDir',fullfile(root,'results','logs',g.outputTag));
cfg.plot.visible=false;cfg.plot.closeAfterExport=true;
cfg.plot.figurePosition=g.plotPosition;cfg.plot.fontSize=g.plotFontSize;
cfg=uav.config.validate(cfg);
indices=g.pointIndices;
assert(all(indices<=numel(archived.cfg.experiments.snapshotSweep)));
points=archived.cfg.experiments.snapshotSweep(indices);
names=["oracle","dl","esp","wco_archived","current_covfit", ...
    "same_nnls0","yue_same_ls","lin_same_literal","lin_same_loaded_noise", ...
    "peak_ridge","peak_nnls0","yue_peak_ls","lin_peak_literal","lin_peak_loaded_noise", ...
    "oracle_incm_nominal","oracle_nodes_nnls0","oracle_nodes_ridge","oracle_noise_current", ...
    "low_load_ridge","low_load_nnls0","low_peak_nnls0","refined_nnls0","yue_refined_ls","lin_refined_literal"];
labels=["理想上界","DL","ESP","WCO（原图）","当前CovFit", ...
    "同候选NNLS零正则","Yue-LS同候选","Lin原式同候选","Lin加载补噪同候选", ...
    "局部峰＋原正则","局部峰＋零正则","Yue-LS局部峰","Lin原式局部峰","Lin加载补噪局部峰", ...
    "诊断：真实INCM＋名义导向","诊断：真候选零正则","诊断：真候选原正则","诊断：真实噪声", ...
    "低加载＋原正则","低加载＋零正则","低加载局部峰","亚网格峰＋NNLS","Yue-LS细化候选","Lin原式细化候选"];
T=cfg.monteCarlo.nTrials;K=cfg.experiments.periodCounts.snapshotSweep;
P=numel(points);J=numel(g.jitterStdDeg);Q=numel(names);
shape=[T,K,P,J,Q];
metricFields={'sinrLinear','desiredGain','interferenceLeakage','incmNmse','weightNorm','weightLoading','negativePowerCount'};
for f=metricFields,raw.(f{1})=nan(shape);end
raw.status=zeros(shape,'int8');raw.errors=strings(shape);
diagnosticFields={'originalMinCoherence','peakMinCoherence','originalCoveredSources', ...
    'peakCoveredSources','originalCandidateCount','peakCandidateCount','noiseEstimate', ...
    'originalPowerRatio','zeroRegularizationPowerRatio','peakPowerRatio','replayAbsoluteError', ...
    'lowLoadMinCoherence','refinedMinCoherence','refinedCandidateCount'};
for f=diagnosticFields,diagnostics.(f{1})=nan(T,K,P,J);end
raw.snapshotSignatures=zeros(T,K,P,J,4);
seeds=cell(P,J);
sourceDesired=cfg.attitude.stdDeg;
for ji=1:J
    pointCfg=cfg;pointCfg.attitude.stdDeg(:)=g.jitterStdDeg(ji);
    pointCfg=uav.config.validate(pointCfg);
    for ip=1:P
        pointCfg.timeline.L=points(ip);
        pointCfg=uav.config.validate(pointCfg);
        % 原图使用实验7、原扫描点序号；稀疏抽点不重编号，保证精确复放。
        originalPoint=indices(ip);
        seed=archived.seedScheme{originalPoint};seeds{ip,ji}=seed;
        for tr=1:T
            scene=uav.scenario.generate(pointCfg,uav.util.makeStream(cfg.random.baseSeed,seed.scenarioSubstreams(tr)));
            for k=1:K
                obs=uav.scenario.generatePeriod(pointCfg,scene,k, ...
                    uav.util.makeStream(cfg.random.baseSeed,seed.snapshotSubstreams(tr,k)));
                safe=uav.experiment.algorithmObservation(obs);
                [originalR,od]=uav.incm.upaCovfit(safe,pointCfg);
                A=od.candidateManifold;D=size(od.desiredCandidateNodes,1);noise=od.noisePower;
                [peakA,peakD]=uav.spectrum.localPeakCandidates(safe.Rx,od,pointCfg);
                lowCfg=pointCfg;lowCfg.incm.staticCovfit.loadingRatio=g.lowLoadingRatio;
                [~,lowD]=uav.incm.upaCovfit(safe,lowCfg);
                lowA=lowD.candidateManifold;lowDesired=size(lowD.desiredCandidateNodes,1);
                [lowPeakA,lowPeakDesired]=uav.spectrum.localPeakCandidates(safe.Rx,lowD,lowCfg);
                [refinedA,~,refinement]=uav.spectrum.refineInterfererCandidates(safe,lowD,lowCfg);
                trueA=obs.truth.trueASVs;
                [diagnostics.originalMinCoherence(tr,k,ip,ji),diagnostics.originalCoveredSources(tr,k,ip,ji)]=coverage(A(:,D+1:end),trueA(:,2:end),g.coverageCoherence);
                [diagnostics.peakMinCoherence(tr,k,ip,ji),diagnostics.peakCoveredSources(tr,k,ip,ji)]=coverage(peakA(:,peakD+1:end),trueA(:,2:end),g.coverageCoherence);
                diagnostics.originalCandidateCount(tr,k,ip,ji)=size(A,2)-D;
                diagnostics.peakCandidateCount(tr,k,ip,ji)=size(peakA,2)-peakD;
                diagnostics.noiseEstimate(tr,k,ip,ji)=noise;
                diagnostics.lowLoadMinCoherence(tr,k,ip,ji)=coverage(lowA(:,lowDesired+1:end),trueA(:,2:end),g.coverageCoherence);
                diagnostics.refinedMinCoherence(tr,k,ip,ji)=coverage(refinedA(:,refinement.desiredCount+1:end),trueA(:,2:end),g.coverageCoherence);
                diagnostics.refinedCandidateCount(tr,k,ip,ji)=size(refinedA,2)-refinement.desiredCount;
                diagnostics.originalPowerRatio(tr,k,ip,ji)=sum(od.interfererFittedPowers)/sum(obs.truth.interfererPowers);
                raw.snapshotSignatures(tr,k,ip,ji,:)=uav.util.snapshotSignature(obs.X);
                for mi=1:Q
                    Rh=[];pd=struct('negativePowerCount',NaN);
                    try
                        switch names(mi)
                            case "oracle"
                                [w,wd]=uav.beamforming.solve('mvdr',obs.truth.trueIncm,trueA(:,1),pointCfg);
                            case "dl"
                                [w,wd]=uav.beamforming.solve('dlmvdr',safe.Rx,cfg.derived.nominalSteeringReference,pointCfg);
                            case "esp"
                                [w,wd]=uav.beamforming.eigenspaceProjection(safe.Rx,cfg.derived.nominalSteeringReference,pointCfg);
                            case "wco_archived"
                                if g.jitterStdDeg(ji)~=sourceDesired(1)
                                    raw.errors(tr,k,ip,ji,mi)="零抖动诊断不复算WCO；原图归档只适用于5°";continue;
                                end
                                ai=find(archived.methodNames=="wco",1);
                                raw.sinrLinear(tr,k,ip,ji,mi)=archived.raw.sinrLinear(tr,k,originalPoint,ai);
                                raw.desiredGain(tr,k,ip,ji,mi)=archived.raw.desiredGain(tr,k,originalPoint,ai);
                                raw.interferenceLeakage(tr,k,ip,ji,mi)=archived.raw.interferenceLeakage(tr,k,originalPoint,ai);
                                raw.status(tr,k,ip,ji,mi)=1;continue;
                            case "current_covfit"
                                Rh=originalR;
                            case "oracle_incm_nominal"
                                Rh=obs.truth.trueIncm;
                            otherwise
                                dictionary=A;desired=D;estimateNoise=noise;rule="nnls";lambda=0;
                                if any(names(mi)==["peak_ridge","peak_nnls0","yue_peak_ls","lin_peak_literal","lin_peak_loaded_noise"])
                                    dictionary=peakA;desired=peakD;
                                elseif any(names(mi)==["oracle_nodes_nnls0","oracle_nodes_ridge"])
                                    dictionary=trueA;desired=1;
                                elseif any(names(mi)==["low_load_ridge","low_load_nnls0"])
                                    dictionary=lowA;desired=lowDesired;
                                elseif names(mi)=="low_peak_nnls0"
                                    dictionary=lowPeakA;desired=lowPeakDesired;
                                elseif any(names(mi)==["refined_nnls0","yue_refined_ls","lin_refined_literal"])
                                    dictionary=refinedA;desired=refinement.desiredCount;
                                elseif names(mi)=="oracle_noise_current"
                                    estimateNoise=cfg.signal.noisePower;
                                end
                                if any(names(mi)==["peak_ridge","oracle_nodes_ridge","oracle_noise_current","low_load_ridge"])
                                    lambda=cfg.incm.staticCovfit.fitRegularization;
                                elseif any(names(mi)==["yue_same_ls","yue_peak_ls","yue_refined_ls"])
                                    rule="yue_ls";
                                elseif any(names(mi)==["lin_same_literal","lin_peak_literal","lin_refined_literal"])
                                    rule="lin_literal";
                                    if size(safe.X,2)<cfg.derived.M
                                        raw.errors(tr,k,ip,ji,mi)="L<M，原式未加载SCM不可逆";continue;
                                    end
                                elseif any(names(mi)==["lin_same_loaded_noise","lin_peak_loaded_noise"])
                                    rule="lin_loaded_noise";
                                end
                                [Rh,pd]=uav.incm.literaturePointPowers(rule,safe.Rx,dictionary,desired,estimateNoise,lambda,pointCfg);
                                if names(mi)=="same_nnls0"
                                    diagnostics.zeroRegularizationPowerRatio(tr,k,ip,ji)=sum(pd.interfererPowers)/sum(obs.truth.interfererPowers);
                                elseif names(mi)=="peak_nnls0"
                                    diagnostics.peakPowerRatio(tr,k,ip,ji)=sum(pd.interfererPowers)/sum(obs.truth.interfererPowers);
                                end
                        end
                        raw.negativePowerCount(tr,k,ip,ji,mi)=pd.negativePowerCount;
                        if ~isempty(Rh)
                            if min(real(eig(Rh))) < -cfg.beamforming.covariancePsdTolerance*max(norm(Rh,2),1)
                                error('uav:study:IndefiniteLs','无约束LS产生非PSD矩阵，未暗中裁剪负功率。');
                            end
                            [w,wd]=uav.beamforming.solve('mvdr',Rh,cfg.derived.nominalSteeringReference,pointCfg);
                        end
                        truth=obs.truth;if ~isempty(Rh),truth.estimatedIncm=Rh;end
                        value=uav.metrics.evaluate(w,truth,pointCfg);
                        raw.sinrLinear(tr,k,ip,ji,mi)=value.sinrLinear;
                        raw.desiredGain(tr,k,ip,ji,mi)=value.desiredGain;
                        raw.interferenceLeakage(tr,k,ip,ji,mi)=value.interferenceLeakage;
                        raw.weightNorm(tr,k,ip,ji,mi)=norm(w);
                        if isfield(wd,'loading'),raw.weightLoading(tr,k,ip,ji,mi)=wd.loading;end
                        raw.negativePowerCount(tr,k,ip,ji,mi)=pd.negativePowerCount;
                        if isfield(value,'incmNmse'),raw.incmNmse(tr,k,ip,ji,mi)=value.incmNmse;end
                        raw.status(tr,k,ip,ji,mi)=1;
                    catch ex
                        raw.status(tr,k,ip,ji,mi)=-1;
                        if strcmp(ex.identifier,'uav:study:NotApplicable'),raw.status(tr,k,ip,ji,mi)=0;end
                        raw.errors(tr,k,ip,ji,mi)=string(ex.identifier)+": "+string(ex.message);
                    end
                end
                active=reshape(raw.sinrLinear(tr,k,ip,ji,:),1,[]);
                assert(all(active(isfinite(active))<=active(1)+cfg.experiments.oracleTolerance),'Oracle上界不应被超过。');
                assert(raw.status(tr,k,ip,ji,5)==1,'当前CovFit复放不能失败。');
                if g.jitterStdDeg(ji)==sourceDesired(1)
                    ai=find(archived.methodNames=="static_incm_mvdr",1);
                    delta=abs(active(5)-archived.raw.sinrLinear(tr,k,originalPoint,ai));
                    diagnostics.replayAbsoluteError(tr,k,ip,ji)=delta;
                    assert(delta<=g.replayTolerance*max(active(5),1),'未精确复放目标图，禁止继续归因。');
                end
            end
        end
        fprintf('Lin/Yue研究：抖动%g°，L=%d，%d/%d完成\n',g.jitterStdDeg(ji),points(ip),ip,P);
    end
end
summaries=cell(1,J);rows=cell(P*J*Q,1);ri=0;
for ji=1:J
    values=reshape(raw.sinrLinear(:,:,:,ji,:),T,K,P,Q);
    summaries{ji}=uav.experiment.summarizePoints(values,cfg,7);
    for ip=1:P
        for mi=1:Q
            ri=ri+1;s=reshape(raw.status(:,:,ip,ji,mi),[],1);
            rows{ri}=table(g.jitterStdDeg(ji),points(ip),names(mi),summaries{ji}.medianDb(ip,mi), ...
                summaries{ji}.ci95Db(1,ip,mi),summaries{ji}.ci95Db(2,ip,mi), ...
                sum(s==1),sum(s==0),sum(s==-1),'VariableNames', ...
                {'jitterStdDeg','snapshotCount','method','medianSinrDb','ciLowDb','ciHighDb','successes','notApplicable','failures'});
        end
    end
end
summaryTable=vertcat(rows{:});
fig=uav.plot.newFigure(cfg);layout=tiledlayout(fig,2,2);
for panel=1:4
    ji=J;if mod(panel,2)==0,ji=1;end
    methods=[1,2,3,4,5,7,8,23,24];
    if panel>2,methods=[1,5,6,10,11,19,20,22,15];end
    if g.jitterStdDeg(ji)==0,methods(methods==4)=[];end
    ax=nexttile(layout);hold(ax,'on');
    for mi=methods
        curve=semilogx(ax,points,summaries{ji}.medianDb(:,mi),'-o','LineWidth',cfg.plot.lineWidth,'MarkerSize',5);
        if mi==1,set(curve,'Color',[0,0,0],'LineStyle','--','Marker','none');end
    end
    grid(ax,'on');xlabel(ax,'每周期快拍数 L');ylabel(ax,'输出SINR中位数（dB）');
    title(ax,sprintf('三轴RMS=%g°',g.jitterStdDeg(ji)));
    legend(ax,cellstr(labels(methods)),'Location','best','FontSize',9);
    set(ax,'FontSize',cfg.plot.fontSize);
end
result=struct('experimentName',"lin_yue_snapshot_study",'methodNames',names,'methodLabels',labels, ...
    'scanPoints',points,'originalPointIndices',indices,'jitterStdDeg',g.jitterStdDeg, ...
    'raw',raw,'diagnostics',diagnostics,'summary',{summaries},'seedScheme',{seeds}, ...
    'sourceArtifact',g.sourceMat,'sourceMetadata',archived.metadata, ...
    'scope',"3×3UPA、固定最终名义导向、共同MVDR；论文重构核心受控移植，非完整算法复现");
counts=struct('trials',T,'periodsPerPoint',K,'points',P,'jitterConditions',J,'methods',Q);
result=uav.experiment.finalize(result,summaryTable,fig,cfg,counts,timer);
end

function [minimum,covered]=coverage(candidates,trueInterferers,threshold)
%COVERAGE 真值仅在研究评价层使用，不返回给候选或功率估计器。
if isempty(trueInterferers),minimum=1;covered=0;return;end
if isempty(candidates),minimum=0;covered=0;return;end
coherence=abs(candidates'*trueInterferers).^2/size(candidates,1)^2;
best=max(coherence,[],1);minimum=min(best);covered=sum(best>=threshold);
end
