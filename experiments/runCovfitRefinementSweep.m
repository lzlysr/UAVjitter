function result = runCovfitRefinementSweep(profile)
%RUNCOVFITREFINEMENTSWEEP 32快拍下细化参数的配对单因素敏感性实验。
%   四个参数轴共用同一批新观测；每个点只改变一个参数，不选最佳场景。
%   固定最终名义期望导向；Le完整版本保留201网格和内部期望候选细化。
arguments
    profile (1,1) string = "smoke"
end
timer=tic;request=config(profile);g=request.experiments.linYueStudy;s=g.refinementSweep;
loaded=load(g.sourceMat,'result');cfg=loaded.result.cfg;
cfg.profile=profile;cfg.experiments.linYueStudy=g;
cfg.timeline.L=s.snapshotCount;cfg.attitude.stdDeg(:)=s.jitterStdDeg;
T=cfg.monteCarlo.nTrials;K=cfg.experiments.periodCounts.snapshotSweep;
if profile=="smoke",T=min(T,g.smokeTrials);K=min(K,g.smokePeriods);
elseif profile~="preliminary",error('uav:study:Profile','仅支持smoke/preliminary。');end
cfg.monteCarlo.nTrials=T;cfg.experiments.periodCounts.snapshotSweep=K;
root=fileparts(fileparts(mfilename('fullpath')));
cfg.output=struct('rootDir',root,'overwrite',true, ...
    'rawDir',fullfile(root,'results','raw',g.outputTag), ...
    'tableDir',fullfile(root,'results','tables',g.outputTag), ...
    'figureDir',fullfile(root,'results','figures',g.outputTag), ...
    'logDir',fullfile(root,'results','logs',g.outputTag));
cfg.plot.visible=false;cfg.plot.closeAfterExport=true;
cfg.plot.figurePosition=g.plotPosition;cfg.plot.fontSize=g.plotFontSize;
cfg=uav.config.validate(cfg);
allNames=["oracle","current_covfit","refined_nnls0","le_refined_covfit", ...
    "le_shared_candidates","yue_refined_ls","lin_refined_literal"];
allLabels=["理想上界","原CovFit","本次局部细化CovFit", ...
    "Le完整细化CovFit（201网格）","Le细化器（共同候选）", ...
    "Yue-LS＋本次细化","Lin点功率＋本次细化"];
[known,labelIndices]=ismember(s.methodNames,allNames);assert(all(known),'扫描方法名非法。');
names=s.methodNames;labels=allLabels(labelIndices);
axesNames=["radiusBins","maxIterations","optimalityTolerance","loadingRatio"];
axisLabels=["局部细化半径（网格间隔数；Le不受此参数影响）", ...
    "细化最大迭代次数","细化最优性容差","细化前端谱加载比例"];
points=cell(1,numel(axesNames));
for ai=1:numel(axesNames),points{ai}=s.(axesNames(ai));end
P=max(cellfun(@numel,points));Q=numel(names);shape=[T,K,P,numel(axesNames),Q];
fields={'sinrLinear','desiredGain','interferenceLeakage','incmNmse', ...
    'relativeKktResidual','candidateCount','acceptedFraction','weightLoading'};
for f=fields,raw.(f{1})=nan(shape);end
raw.status=zeros(shape,'int8');raw.errors=strings(shape);
% 观测只生成一次：参数扫描不改变姿态、训练序列、干扰或噪声抽样。
observations=cell(T,K);seeds=struct('scenarioSubstreams',zeros(T,1),'snapshotSubstreams',zeros(T,K));
signatures=zeros(T,K,4);
for tr=1:T
    seeds.scenarioSubstreams(tr)=uav.experiment.substream(cfg,s.experimentNumber,s.seedPointIndex,tr,0,'scenario');
    scene=uav.scenario.generate(cfg,uav.util.makeStream(cfg.random.baseSeed,seeds.scenarioSubstreams(tr)));
    for k=1:K
        seeds.snapshotSubstreams(tr,k)=uav.experiment.substream(cfg,s.experimentNumber,s.seedPointIndex,tr,k,'snapshot');
        observations{tr,k}=uav.scenario.generatePeriod(cfg,scene,k, ...
            uav.util.makeStream(cfg.random.baseSeed,seeds.snapshotSubstreams(tr,k)));
        signatures(tr,k,:)=uav.util.snapshotSignature(observations{tr,k}.X);
    end
end
summaries=cell(1,numel(axesNames));tables={};
for ai=1:numel(axesNames)
    for ip=1:numel(points{ai})
        c=cfg;value=points{ai}(ip);
        switch axesNames(ai)
            case "radiusBins",c.experiments.linYueStudy.refineCandidateRadiusBins=value;
            case "maxIterations"
                c.experiments.linYueStudy.localRefinement.maxIterations=value;
                c.experiments.linYueStudy.leRefinement.maxIterations=value;
            case "optimalityTolerance"
                c.experiments.linYueStudy.localRefinement.optimalityTolerance=value;
                c.experiments.linYueStudy.leRefinement.optimalityTolerance=value;
            case "loadingRatio"
                c.experiments.linYueStudy.lowLoadingRatio=value;
                c.experiments.linYueStudy.leCovfit.loadingRatio=value;
        end
        low=c;low.incm.staticCovfit.loadingRatio=c.experiments.linYueStudy.lowLoadingRatio;
        le=c;le.spectrum.gridSize=g.leGridSize;le=uav.config.validate(le);
        for tr=1:T
            for k=1:K
                obs=observations{tr,k};safe=uav.experiment.algorithmObservation(obs);
                [originalR,original]=uav.incm.upaCovfit(safe,c);
                [~,lowD]=uav.incm.upaCovfit(safe,low);
                [fineA,~,fine]=uav.spectrum.refineInterfererCandidates(safe,lowD,low);
                for mi=1:Q
                    R=[];pd=struct;A=[];lambda=0;noise=lowD.noisePower;accepted=NaN;
                    try
                        switch names(mi)
                            case "oracle"
                                R=obs.truth.trueIncm;
                                [w,wd]=uav.beamforming.solve('mvdr',R,obs.truth.trueASVs(:,1),c);
                            case "current_covfit"
                                R=originalR;A=original.candidateManifold;
                                pd.powers=original.fittedPowers;noise=original.noisePower;
                                lambda=c.incm.staticCovfit.fitRegularization;
                            case "le_refined_covfit"
                                [R,d]=uav.incm.upaCovfitLe(safe,le);
                                A=d.candidateManifold;pd.powers=d.fittedPowers;noise=d.noisePower;
                                flags=[d.candidateRefinement.desired.accepted;d.candidateRefinement.interferer.accepted];
                                accepted=mean(flags);
                            otherwise
                                A=fineA;D=fine.desiredCount;rule="nnls";accepted=mean(fine.accepted);
                                if names(mi)=="le_shared_candidates"
                                    [A,D,d]=uav.spectrum.leSharedCandidates(safe,lowD,low);accepted=mean(d.accepted);
                                elseif names(mi)=="yue_refined_ls",rule="yue_ls";
                                elseif names(mi)=="lin_refined_literal",rule="lin_literal";end
                                [R,pd]=uav.incm.literaturePointPowers(rule,safe.Rx,A,D,noise,lambda,c);
                        end
                        if any(names(mi)==["current_covfit","refined_nnls0","le_refined_covfit","le_shared_candidates"])
                            residual=uav.incm.nonnegativeFitKkt(safe.Rx,A,noise,pd.powers,lambda);
                            raw.relativeKktResidual(tr,k,ip,ai,mi)=residual;
                            if residual>g.auditKktTolerance,error('uav:study:InvalidNnls','非负拟合返回点未通过KKT检查。');end
                        end
                        if names(mi)~="oracle"
                            if min(real(eig(R))) < -c.beamforming.covariancePsdTolerance*max(norm(R,2),1)
                                error('uav:study:IndefiniteLs','重构矩阵非PSD，保留失败而不裁剪功率。');
                            end
                            [w,wd]=uav.beamforming.solve('mvdr',R,c.derived.nominalSteeringReference,c);
                        end
                        truth=obs.truth;truth.estimatedIncm=R;metric=uav.metrics.evaluate(w,truth,c);
                        for f={'sinrLinear','desiredGain','interferenceLeakage','incmNmse'}
                            raw.(f{1})(tr,k,ip,ai,mi)=metric.(f{1});
                        end
                        raw.candidateCount(tr,k,ip,ai,mi)=size(A,2);
                        raw.acceptedFraction(tr,k,ip,ai,mi)=accepted;
                        raw.weightLoading(tr,k,ip,ai,mi)=wd.loading;
                        raw.status(tr,k,ip,ai,mi)=1;
                    catch ex
                        raw.status(tr,k,ip,ai,mi)=-1;
                        raw.errors(tr,k,ip,ai,mi)=string(ex.identifier)+": "+string(ex.message);
                    end
                end
                oracle=find(names=="oracle",1);
                if ~isempty(oracle)
                    values=reshape(raw.sinrLinear(tr,k,ip,ai,:),1,[]);
                    assert(all(values(isfinite(values))<=values(oracle)+c.experiments.oracleTolerance),'超过Oracle。');
                end
            end
        end
        fprintf('32快拍细化扫描：%s=%g，%d/%d完成\n',axesNames(ai),value,ip,numel(points{ai}));
    end
    count=numel(points{ai});values=reshape(raw.sinrLinear(:,:,1:count,ai,:),T,K,count,Q);
    summaries{ai}=uav.experiment.summarizePoints(values,cfg,s.experimentNumber);
    for ip=1:count
        for mi=1:Q
            st=raw.status(:,:,ip,ai,mi);summary=summaries{ai};
            tables{end+1}=table(axesNames(ai),points{ai}(ip),names(mi), ...
                summary.medianDb(ip,mi),summary.ci95Db(1,ip,mi),summary.ci95Db(2,ip,mi), ...
                nnz(st==1),nnz(st==-1),'VariableNames',{'parameter','value','method', ...
                'medianSinrDb','ciLowDb','ciHighDb','successes','failures'}); %#ok<AGROW>
        end
    end
end
result=struct('experimentName',"covfit_refinement_parameter_sweep",'methodNames',names, ...
    'methodLabels',labels,'axisNames',axesNames,'axisLabels',axisLabels,'scanPoints',{points}, ...
    'summary',{summaries},'raw',raw,'seedScheme',seeds,'snapshotSignatures',signatures, ...
    'sourceArtifact',g.sourceMat,'sourceMetadata',loaded.result.metadata, ...
    'scope',"32快拍、5°RMS、逐轴单因素、同观测、固定最终名义期望导向；不据此选择最优参数");
result.cfg=cfg;fig=plotCovfitRefinementSweep(result);
counts=struct('trials',T,'periods',K,'parameterAxes',numel(axesNames),'parameterPoints',sum(cellfun(@numel,points)),'methods',Q);
result=uav.experiment.finalize(result,vertcat(tables{:}),fig,cfg,counts,timer);
end
