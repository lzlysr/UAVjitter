function result=runUpaIntegralGridConvergence(cfg)
%RUNUPAINTEGRALGRIDCONVERGENCE 回放实验8，隔离原始Rec的二维数值求积误差。
%   本诊断不改变正式实验，也不运行processed方法或调整任何算法参数。
%   sample：回放原始有限快拍SCM；population：同一姿态下用理论总协方差；
%   nominal_population：零抖动理论总协方差，用于检查静态积分自身的偏差。
%   后两类明确属于真值诊断，不作为可部署方法，不向sample重构泄漏真值。
timer=tic; settings=cfg.experiments.upaGridConvergence;
sourcePath=fullfile(cfg.output.rootDir,settings.sourceResult);
loaded=load(sourcePath,'result');source=loaded.result;c=source.cfg;
if c.array.Nx<=1 || c.array.Ny<=1 || source.scanKind~="snr"
    error('uav:experiment:UpaSourceRequired','诊断输入必须为UPA的实验8原始结果。');
end
% 物理参数以待解释结果的完整配置为准；本次新增诊断设置来自当前config。
c=uav.config.validate(c);c.experiments.upaGridConvergence=settings;
trials=settings.trialIndices; periods=settings.periodIndices;
if max(trials)>size(source.raw.sinrLinear,1) || max(periods)>size(source.raw.sinrLinear,2)
    error('uav:experiment:ReplayRange','诊断试次/周期超出已保存实验的范围，请在config设置。');
end
[found,pointIndices]=ismember(settings.snrDb,source.scanPoints);
if ~all(found),error('uav:experiment:ReplaySnr','诊断SNR不在原实验扫描点中。');end
M=c.derived.M; count=numel(trials)*numel(periods)*numel(pointIndices);
total=2*count+numel(pointIndices);
covariances=complex(zeros(M,M,total));truths=cell(total,1);
kind=strings(total,1);snr=zeros(total,1);trial=zeros(total,1);period=zeros(total,1);
oracle=zeros(total,1);fixedTrue=zeros(total,1); savedOrigin=nan(total,1);
savedPlain=nan(total,1);savedRefined=nan(total,1);signatureError=zeros(total,1);
originIndex=find(source.methodNames=="incm_integral_origin",1);
plainIndex=find(source.methodNames=="incm_covfit",1);
refinedIndex=find(source.methodNames=="incm_covfit_refined",1);
b=0;
for p=pointIndices
    pc=c;pc.signal.snrDb=source.scanPoints(p);pc=uav.config.validate(pc);
    for t=trials
        sceneSeed=uav.experiment.substream(pc,8,p,t,0,"scenario");
        scene=uav.scenario.generate(pc,uav.util.makeStream(pc.random.baseSeed,sceneSeed));
        for k=periods
            seed=uav.experiment.substream(pc,8,p,t,k,"snapshot");
            obs=uav.scenario.generatePeriod(pc,scene,k,uav.util.makeStream(pc.random.baseSeed,seed));
            b=b+1;populate(b,"sample",obs.Rx,obs.truth,pc,t,k);
            actual=uav.util.snapshotSignature(obs.X);
            expected=reshape(source.raw.snapshotSignatures(t,k,p,originIndex,:),1,[]);
            signatureError(b)=max(abs(actual-expected));
            savedOrigin(b)=source.raw.sinrLinear(t,k,p,originIndex);
            if ~isempty(plainIndex),savedPlain(b)=source.raw.sinrLinear(t,k,p,plainIndex);end
            if ~isempty(refinedIndex),savedRefined(b)=source.raw.sinrLinear(t,k,p,refinedIndex);end
            population=obs.truth.trueIncm+obs.truth.desiredPower* ...
                obs.truth.trueASVs(:,1)*obs.truth.trueASVs(:,1)';
            populate(b+count,"population",population,obs.truth,pc,t,k);
        end
    end
    zero=pc;zero.attitude.stdDeg=[0,0,0];zero=uav.config.validate(zero);
    scene=uav.scenario.generate(zero,uav.util.makeStream(zero.random.baseSeed,sceneSeed));
    obs=uav.scenario.generatePeriod(zero,scene,1,uav.util.makeStream(zero.random.baseSeed,seed));
    population=obs.truth.trueIncm+obs.truth.desiredPower*obs.truth.trueASVs(:,1)*obs.truth.trueASVs(:,1)';
    index=2*count+find(pointIndices==p,1);
    populate(index,"nominal_population",population,obs.truth,zero,0,0);
end
if any(signatureError~=0)
    error('uav:experiment:ReplayMismatch','快拍签名与已保存实验不一致，禁止继续比较。');
end
gridSizes=[settings.gridSizes,settings.midpointGridSizes];
rules=[repmat("endpoint",size(settings.gridSizes)),repmat("midpoint",size(settings.midpointGridSizes))];
G=numel(gridSizes);Rall=complex(zeros(M,M,total,G));
sinrLinear=zeros(total,G);incmNmse=zeros(total,G);shapeError=zeros(total,G);
leakage=zeros(total,G);desiredGain=zeros(total,G);loading=zeros(total,G);
quadrature=cell(G,1);gridRuntime=zeros(G,1);
previousThreads=maxNumCompThreads(settings.threadCount);
restore=onCleanup(@()maxNumCompThreads(previousThreads));
cached=[];reused=false(G,1);
if isfield(settings,'reuseResult') && strlength(settings.reuseResult)>0
    old=load(settings.reuseResult,'result');cached=old.result;
    % 只有同一组协方差、同一阵列/积分域时才允许复用；不按文件名猜测一致。
    if ~isequal(cached.inputCovariances,covariances) || ~isequal(cached.cfg.array,c.array) || ...
            ~isequal(cached.cfg.derived.nominalSteeringReference,c.derived.nominalSteeringReference) || ...
            cached.cfg.spectrum.desiredExclusionHalfWidthDeg~=c.spectrum.desiredExclusionHalfWidthDeg
        error('uav:experiment:InvalidGridCache','旧诊断输入或积分域不一致，不允许复用。');
    end
end
for g=1:G
    gridTimer=tic;
    match=[];
    if ~isempty(cached),match=find(cached.gridSizes==gridSizes(g) & cached.rules==rules(g),1);end
    if isempty(match)
        [Rall(:,:,:,g),quadrature{g}]=uav.incm.integrateUpaBatch(covariances,c,gridSizes(g),rules(g),settings);
    else
        Rall(:,:,:,g)=cached.reconstructedIncm(:,:,:,match);
        quadrature{g}=cached.quadrature{match};reused(g)=true;
    end
    for n=1:total
        R=Rall(:,:,n,g);
        [w,d]=uav.beamforming.solve("mvdr",R,c.derived.nominalSteeringReference,c);
        truth=truths{n};truth.estimatedIncm=R;
        m=uav.metrics.evaluate(w,truth,c);
        if m.sinrLinear>oracle(n)*(1+c.experiments.oracleTolerance)
            error('uav:experiment:OracleBoundViolation','诊断SINR超过真实最优上界。');
        end
        sinrLinear(n,g)=m.sinrLinear;incmNmse(n,g)=m.incmNmse;
        desiredGain(n,g)=m.desiredGain;leakage(n,g)=m.interferenceLeakage;loading(n,g)=d.loading;
        T=truth.trueIncm;
        shapeError(n,g)=norm(R/real(trace(R))-T/real(trace(T)),'fro');
    end
    gridRuntime(g)=toc(gridTimer);
    if reused(g),gridRuntime(g)=cached.summary.gridRuntimeSec(find(cached.summary.gridSize==gridSizes(g) & cached.summary.rule==rules(g),1));end
    fprintf('UPA原始Rec积分收敛：%s G=%d，节点%d，耗时%.2f秒\n', ...
        rules(g),gridSizes(g),quadrature{g}.includedNodeCount,gridRuntime(g));
end
sourceGrid=find(gridSizes==c.spectrum.gridSize & rules=="endpoint",1);
if isempty(sourceGrid),error('uav:experiment:MissingSourceGrid','扫描必须包含原结果网格。');end
selected=kind=="sample";
replayError=max(abs(10*log10(sinrLinear(selected,sourceGrid)./savedOrigin(selected))));
if replayError>settings.replayToleranceDb
    error('uav:experiment:ReplaySinrMismatch','原始网格SINR回放不一致：%.3g dB。',replayError);
end
% 独立调用生产实现校验代数加速，防止“加速后的另一个算法”误导诊断。
grid=c.derived.spectrumGrid;keep=uav.spectrum.admissibleDirections(grid.directionCosines.',c);
direct=uav.incm.guLeshem2012Rec(covariances(:,:,1),grid.manifold(:,keep),grid.weights(keep));
parityError=norm(direct-Rall(:,:,1,sourceGrid),'fro')/norm(direct,'fro');
if parityError>settings.parityTolerance,error('uav:experiment:IntegralParity','加速积分与生产实现不一致。');end
reference=find(rules=="endpoint",1,'last');
rows=cell(0,21);cases=["sample","population","nominal_population"];
for caseName=cases
    for inputSnr=settings.snrDb
        indices=find(kind==caseName & snr==inputSnr);
        for g=1:G
            x=sinrLinear(indices,g);delta=abs(10*log10(x./sinrLinear(indices,reference)));
            rel=zeros(numel(indices),1);prev=nan(numel(indices),1);prevDb=prev;
            for j=1:numel(indices)
                n=indices(j);R=Rall(:,:,n,g);target=Rall(:,:,n,reference);
                rel(j)=norm(R-target,'fro')/norm(target,'fro');
                if g>1 && rules(g)=="endpoint"
                    last=Rall(:,:,n,g-1);prev(j)=norm(R-last,'fro')/norm(R,'fro');
                    prevDb(j)=abs(10*log10(sinrLinear(n,g)/sinrLinear(n,g-1)));
                end
            end
            rows(end+1,:)={caseName,inputSnr,rules(g),gridSizes(g),quadrature{g}.step, ...
                quadrature{g}.includedNodeCount,quadrature{g}.includedArea,numel(indices), ...
                10*log10(median(x)),10*log10(quantile(x,.05)),10*log10(median(oracle(indices))), ...
                10*log10(median(fixedTrue(indices))),median(incmNmse(indices,g)), ...
                median(shapeError(indices,g)),median(leakage(indices,g)),max(delta), ...
                quantile(delta,.95),max(rel),max(prevDb),max(prev),gridRuntime(g)}; %#ok<AGROW>
        end
    end
end
summary=cell2table(rows,'VariableNames',{'caseName','inputSnrDb','rule','gridSize','step', ...
    'includedNodes','includedArea','sampleCount','medianSinrDb','p05SinrDb','oracleMedianDb', ...
    'trueIncmFixedReferenceMedianDb','medianIncmNmse','medianTraceNormalizedShapeError', ...
    'medianInterferenceLeakage','maxAbsSinrDeltaToFinestDb','p95AbsSinrDeltaToFinestDb', ...
    'maxRelativeMatrixDeltaToFinest','maxAbsSinrDeltaToPreviousDb','maxRelativeMatrixDeltaToPrevious','gridRuntimeSec'});
% 最后两次端点加密及最高分辨率半格错位分别检查；不把参考网格自身的零误差当收敛证据。
checkIndices=[reference-1,reference,find(rules=="midpoint",1,'last')];
checks=summary(ismember(summary.gridSize,gridSizes(checkIndices)),:);
passed=true(height(checks),1);
for j=1:height(checks)
    if checks.rule(j)=="endpoint"
        passed(j)=checks.maxAbsSinrDeltaToPreviousDb(j)<=settings.sinrToleranceDb && ...
            checks.maxRelativeMatrixDeltaToPrevious(j)<=settings.matrixRelativeTolerance;
    else
        passed(j)=checks.maxAbsSinrDeltaToFinestDb(j)<=settings.sinrToleranceDb && ...
            checks.maxRelativeMatrixDeltaToFinest(j)<=settings.matrixRelativeTolerance;
    end
end
checks.passed=passed;
result=struct('cfg',c,'settings',settings,'sourceMetadata',source.metadata, ...
    'sourcePath',sourcePath,'sampleIndex',table(kind,snr,trial,period), ...
    'inputCovariances',covariances,'truth',{truths},'gridSizes',gridSizes,'rules',rules, ...
    'quadrature',{quadrature},'reconstructedIncm',Rall,'sinrLinear',sinrLinear, ...
    'incmNmse',incmNmse,'traceNormalizedShapeError',shapeError, ...
    'desiredGain',desiredGain,'interferenceLeakage',leakage,'mvdrNumericalLoading',loading, ...
    'oracleLinear',oracle,'trueIncmFixedReferenceLinear',fixedTrue, ...
    'savedOriginLinear',savedOrigin,'savedCovfitLinear',savedPlain,'savedRefinedLinear',savedRefined, ...
    'summary',summary,'convergenceChecks',checks,'allChecksPassed',all(passed), ...
    'reusedGridResults',reused, ...
    'replayMaxErrorDb',replayError,'productionParityRelativeError',parityError);
result.metadata=uav.experiment.metadata(c,"diagnosis_upa_integral_grid_convergence", ...
    struct('independentTrials',numel(trials),'periodsPerTrial',numel(periods), ...
    'snrPoints',numel(pointIndices),'sampleFrames',count,'totalCovariances',total),toc(timer));
stamp=char(datetime('now','TimeZone','Asia/Shanghai','Format','yyyyMMdd_HHmmss'));
name=['diagnosis_upa_integral_grid_' stamp];
result.artifacts=struct('mat',fullfile(c.output.rawDir,[name '.mat']), ...
    'summaryCsv',fullfile(c.output.tableDir,[name '_summary.csv']), ...
    'checksCsv',fullfile(c.output.tableDir,[name '_checks.csv']));
save(result.artifacts.mat,'result','-v7.3');
writetable(summary,result.artifacts.summaryCsv);writetable(checks,result.artifacts.checksCsv);
fprintf('回放最大误差=%.3g dB，生产积分相对误差=%.3g，全部收敛检查=%d\n',replayError,parityError,all(passed));
disp(result.artifacts);
    function populate(index,caseName,R,truth,pc,t,k)
        covariances(:,:,index)=R;truths{index}=truth;kind(index)=caseName;
        snr(index)=pc.signal.snrDb;trial(index)=t;period(index)=k;
        w=uav.beamforming.solve("mvdr",truth.trueIncm,truth.trueASVs(:,1),pc);
        metric=uav.metrics.evaluate(w,truth,pc);oracle(index)=metric.sinrLinear;
        w=uav.beamforming.solve("mvdr",truth.trueIncm,pc.derived.nominalSteeringReference,pc);
        metric=uav.metrics.evaluate(w,truth,pc);fixedTrue(index)=metric.sinrLinear;
    end
end
