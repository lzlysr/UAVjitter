function audit = auditCovfitRefinementDifferences(matPath)
%AUDITCOVFITREFINEMENTDIFFERENCES 复放全部大于预定0.1dB差异的样本作事后归因。
%   不删除样本、不修改性能主表；真值仅评价候选覆盖和固定先验分类。
%   审计集合由差异阈值确定，不能用它估计整个场景的无条件性能。
loaded=load(matPath,'result');r=loaded.result;cfg=r.cfg;g=cfg.experiments.linYueStudy;
localIndex=find(r.methodNames=="refined_nnls0");leIndex=find(r.methodNames=="le_refined_covfit");
a=r.raw.sinrLinear(:,:,:,:,localIndex);b=r.raw.sinrLinear(:,:,:,:,leIndex);
delta=10*log10(b./a);selected=find(abs(delta)>g.comparisonToleranceDb);
rows=cell(numel(selected),1);
for q=1:numel(selected)
    [tr,k,ip,ji]=ind2sub(size(delta),selected(q));
    c=cfg;c.timeline.L=r.scanPoints(ip);c.attitude.stdDeg(:)=r.jitterStdDeg(ji);c=uav.config.validate(c);
    seed=r.seedScheme{ip,ji};
    scene=uav.scenario.generate(c,uav.util.makeStream(c.random.baseSeed,seed.scenarioSubstreams(tr)));
    obs=uav.scenario.generatePeriod(c,scene,k,uav.util.makeStream(c.random.baseSeed,seed.snapshotSubstreams(tr,k)));
    safe=uav.experiment.algorithmObservation(obs);
    low=c;low.incm.staticCovfit.loadingRatio=g.lowLoadingRatio;
    [~,d]=uav.incm.upaCovfit(safe,low);
    [A,~,fine]=uav.spectrum.refineInterfererCandidates(safe,d,low);
    [R,power]=uav.incm.literaturePointPowers('nnls',safe.Rx,A,fine.desiredCount,d.noisePower,0,c);
    le=c;le.spectrum.gridSize=g.leGridSize;le=uav.config.validate(le);
    [S,legacy]=uav.incm.upaCovfitLe(safe,le);
    [w,~]=uav.beamforming.solve('mvdr',R,c.derived.nominalSteeringReference,c);
    [v,~]=uav.beamforming.solve('mvdr',S,c.derived.nominalSteeringReference,c);
    ma=uav.metrics.evaluate(w,obs.truth,c);mb=uav.metrics.evaluate(v,obs.truth,c);
    replay=max(abs([ma.sinrLinear-a(selected(q)),mb.sinrLinear-b(selected(q))]));
    assert(replay<=g.replayTolerance*max([a(selected(q)),b(selected(q)),1]),'差异样本未准确复放。');
    outside=uav.spectrum.admissibleDirections(obs.truth.directionCosines.',c);
    [ca,na]=coverage(A(:,fine.desiredCount+1:end),obs.truth.trueASVs(:,2:end),g.coverageCoherence);
    D=size(legacy.desiredCandidateNodes,1);
    [cb,nb]=coverage(legacy.candidateManifold(:,D+1:end),obs.truth.trueASVs(:,2:end),g.coverageCoherence);
    ka=uav.incm.nonnegativeFitKkt(safe.Rx,A,d.noisePower,power.powers,0);
    kb=uav.incm.nonnegativeFitKkt(safe.Rx,legacy.candidateManifold,legacy.noisePower,legacy.fittedPowers,0);
    % 额外诊断：固定Le干扰节点，只将拟合期望字典限制为名义单点。
    % 不使用真实方向；不是新增主方法，也不替换原性能结果。
    controlNodes=[c.derived.nominalDirectionCosines(:,1).';legacy.interfererCandidateNodes];
    controlA=uav.array.steeringUPA(controlNodes.',c.array.Nx,c.array.Ny,c.array.spacingOverLambda);
    controlR=uav.incm.literaturePointPowers('nnls',safe.Rx,controlA,1,legacy.noisePower,0,c);
    controlW=uav.beamforming.solve('mvdr',controlR,c.derived.nominalSteeringReference,c);
    control=uav.metrics.evaluate(controlW,obs.truth,c);
    trueI=obs.truth.trueASVs(:,2:end);
    localProjection=real(sum(conj(trueI).*((R-d.noisePower*eye(c.derived.M))*trueI),1)).'/c.derived.M^2;
    leProjection=real(sum(conj(trueI).*((S-legacy.noisePower*eye(c.derived.M))*trueI),1)).'/c.derived.M^2;
    rows{q}=table(tr,k,c.timeline.L,r.jitterStdDeg(ji),delta(selected(q)),ma.sinrDb,mb.sinrDb, ...
        ~outside(1),nnz(~outside(2:end)),size(A,2)-fine.desiredCount, ...
        legacy.interfererCandidateCount,ca,cb,na,nb,ka,kb, ...
        sum(power.powers(1:fine.desiredCount)),sum(legacy.desiredFittedPowers), ...
        min(localProjection./obs.truth.interfererPowers(:)), ...
        min(leProjection./obs.truth.interfererPowers(:)),control.sinrDb,replay, ...
        'VariableNames',{'trial','period','snapshotCount','jitterStdDeg','leMinusLocalDb', ...
        'localSinrDb','leSinrDb','desiredInsidePrior','interferersInsideDesiredPrior', ...
        'localCandidateCount','leCandidateCount','localMinimumCoherence','leMinimumCoherence', ...
        'localCoveredInterferers','leCoveredInterferers','localKkt','leKkt', ...
        'localDesiredFittedPower','leDesiredFittedPower','localMinInterfererProjectionRatio', ...
        'leMinInterfererProjectionRatio','leNominalOnlyControlSinrDb','replayLinearError'});
end
records=vertcat(rows{:});audit=struct('records',records,'sourceMat',matPath, ...
    'selection',"全部逐周期绝对差异大于配置0.1dB的样本；仅事后诊断", ...
    'totalComparedSamples',numel(delta),'sourceCode',fileread(mfilename('fullpath')+".m"));
prefix="lin_yue_"+cfg.profile+"_refinement_difference_audit";
writetable(records,fullfile(cfg.output.tableDir,prefix+".csv"));
save(fullfile(cfg.output.rawDir,prefix+".mat"),'audit');
fprintf('差异审计完成：%d/%d个样本；期望区包含真实干扰%d次；最大KKT=%g\n', ...
    height(records),numel(delta),nnz(records.interferersInsideDesiredPrior>0), ...
    max([records.localKkt;records.leKkt]));
end

function [minimum,count]=coverage(A,trueA,threshold)
%COVERAGE 只读真值评价，结果不返回候选选择器或拟合器。
if isempty(A),minimum=0;count=0;return;end
best=max(abs(A'*trueA).^2/size(A,1)^2,[],1);
minimum=min(best);count=nnz(best>=threshold);
end
