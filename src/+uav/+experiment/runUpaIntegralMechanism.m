function result=runUpaIntegralMechanism(cfg)
%RUNUPAINTEGRALMECHANISM 用总体协方差区分积分离散误差与估计模型误差。
%   固定上一轮的名义3×3几何，不校正期望导向，不挑选有利随机实现。
%   真实INCM输入、真实干扰点重构、已知噪声替换均明确是Oracle诊断，
%   不能把这些对照冒充不依赖真值的原始Rec或可部署改进方法。
timer=tic;s=cfg.experiments.upaIntegralMechanism;
loaded=load(fullfile(cfg.output.rootDir,s.sourceJacobianResult),'result');old=loaded.result;
c=old.cfg;c.experiments.upaIntegralMechanism=s;
if c.array.Nx<=1 || c.array.Ny<=1 || ~old.allChecksPassed
    error('uav:experiment:InvalidMechanismSource','必须读取已收敛的UPA Jacobian实验。');
end
kernel=cfg.experiments.upaGridConvergence;
previousThreads=maxNumCompThreads(kernel.threadCount);
restore=onCleanup(@()maxNumCompThreads(previousThreads));
ix=zeros(1,numel(s.snrDb));
for k=1:numel(ix)
    ix(k)=find(old.sampleIndex.kind=="nominal_population" & old.sampleIndex.snr==s.snrDb(k),1);
end
M=c.derived.M;B=numel(ix);truth=old.truth(ix);Rx=old.inputCovariances(:,:,ix);
Rin=zeros(M,M,B);
for k=1:B,Rin(:,:,k)=truth{k}.trueIncm;end
noise=c.derived.noisePower*eye(M);
X=cat(3,Rx,Rin(:,:,1),noise);
measures=["dxi_x_dxi_y","solid_angle","angle_product"];
stamp=char(datetime('now','TimeZone','Asia/Shanghai','Format','yyyyMMdd_HHmmss'));
name=['diagnosis_upa_integral_mechanism_' stamp];
checkpointPath=fullfile(c.output.rawDir,[name '_checkpoint.mat']);
records={};labels=cell(0,3);rows=cell(0,12);
% 第一组：换坐标而不换测度，应该收敛到上一轮相同的积分。
for G=s.gridSizes
    [R,d]=uav.incm.integrateUpaBatch(X,c,G,"azel_midpoint",kernel,measures);
    store(R,d,"azel");
end
lastAzel=R;
% 第二组：12801网格只用于上述五个无随机误差的总体场景及控制输入。
% 不把少量场景的极细网格验证冒充全部905个输入均已计算到12801。
[extended,extendedDiagnostics]=uav.incm.integrateUpaBatch(X,c,s.extendedGridSize,"surface_midpoint",kernel,measures(1:2));
store(extended,extendedDiagnostics,"surface_extended");
% 第三组：只改变球冠排除半径。真实INCM输入不包含期望信号，
% 因而其误差不能归因于SCM中的期望信号或有限快拍。
for width=s.exclusionHalfWidthsDeg
    control=c;control.spectrum.desiredExclusionHalfWidthDeg=width;
    [R,d]=uav.incm.integrateUpaBatch(X,control,s.sectorGridSize,"surface_midpoint",kernel,measures(1:2));
    store(R,d,"sector",width);
end
% 独立加密最小的非零排除区，不用默认25度的收敛结论替代它。
control=c;control.spectrum.desiredExclusionHalfWidthDeg=s.sectorVerificationHalfWidthDeg;
[sectorFine,d]=uav.incm.integrateUpaBatch(X,control,max(s.gridSizes),"surface_midpoint",kernel,measures(1:2));
store(sectorFine,d,"sector_verification",s.sectorVerificationHalfWidthDeg);
% 第四组：单个固定干扰，直接比较Capon峰高与积分后的峰面积。
% 干扰方向取当前第一邻区；不把真实方向提供给正式算法。
a=truth{1}.trueASVs(:,2);P=10.^(s.singleSourceInrDb/10)*c.derived.noisePower;
single=zeros(M,M,numel(P)+1);single(:,:,end)=noise;
for k=1:numel(P),single(:,:,k)=P(k)*(a*a')+noise;end
control=c;control.spectrum.desiredExclusionHalfWidthDeg=0;
[singleRec,singleQuadrature]=uav.incm.integrateUpaBatch(single,control,max(s.gridSizes), ...
    "surface_midpoint",kernel,measures(1:2));
singleCoarse=uav.incm.integrateUpaBatch(single,control,s.sectorGridSize, ...
    "surface_midpoint",kernel,measures(1:2));
powerRows=cell(0,8);
for k=1:numel(P)
    pointPower=1/real(a'*(single(:,:,k)\a));
    for m=1:2
        excessMass=real(trace(singleRec(:,:,k,m)-singleRec(:,:,end,m)))/M;
        coarseMass=real(trace(singleCoarse(:,:,k,m)-singleCoarse(:,:,end,m)))/M;
        powerRows(end+1,:)={s.singleSourceInrDb(k),measures(m),P(k),pointPower, ...
            excessMass,excessMass/P(k),singleQuadrature(m).includedArea,abs(coarseMass/excessMass-1)}; %#ok<AGROW>
    end
end
powerTable=cell2table(powerRows,'VariableNames',{'inputInrDb','measure','truePower', ...
    'caponAtTrueDirection','integratedExcessMass','massToTruePower','domainArea','relativeMassGridDelta'});
% 第五组：不进行方向估计，只用真实干扰点读Capon功率的Oracle定位对照。
% 同时做已知噪声底替换，检验噪声各向异性是否能独立解释全部损失。
ablationRows=cell(0,6);
for k=1:B
    Ai=truth{k}.trueASVs(:,2:end);sourcePower=1./real(sum(conj(Ai).*(Rx(:,:,k)\Ai),1));
    point=Ai*diag(sourcePower)*Ai'+noise;
    appendAblation("oracle_point_capon",k,"point",point);
    for m=1:2
        Rraw=extended(:,:,k,m);Rwhite=extended(:,:,end,m);
        % M/域面积仅由白输入的迹归一化确定，不用真实干扰功率调系数。
        calibrated=M/extendedDiagnostics(m).includedArea*(Rraw-Rwhite)+noise;
        appendAblation("oracle_white_baseline_replaced",k,measures(m),calibrated);
    end
    appendAblation("oracle_true_incm",k,"none",truth{k}.trueIncm);
end
ablations=cell2table(ablationRows,'VariableNames',{'controlName','inputSnrDb','measure', ...
    'sinrDb','interferenceLeakage','noiseOutputPower'});
% 仅为根因定位而扫描谱超额部分/白噪声的相对标度；用真值选优的结果
% 是Oracle诊断，不能作为正式算法，也不能据此修改默认配置。
contrastRows=cell(0,6);snrIndex=find(s.snrDb==10,1);
for m=1:2
    for multiplier=s.contrastMultipliers
        R=multiplier*(extended(:,:,snrIndex,m)-extended(:,:,end,m))+noise;
        value=metric(R,snrIndex);
        contrastRows(end+1,:)={measures(m),multiplier,value.sinrDb,value.interferenceLeakage, ...
            value.noiseOutputPower,min(real(eig(R)))}; %#ok<AGROW>
    end
end
contrast=cell2table(contrastRows,'VariableNames',{'measure','multiplier','sinrDb', ...
    'interferenceLeakage','noiseOutputPower','minimumEigenvalue'});
summary=cell2table(rows,'VariableNames',{'inputName','inputSnrDb','groupName','gridSize', ...
    'exclusionHalfWidthDeg','measure','sinrDb','interferenceLeakage','noiseOutputPower', ...
    'topEigenvalue','smallEigenvalueMean','traceNormalizedShapeError'});
reference=find(old.rules=="surface_midpoint",1,'last');checks=cell(0,5);
for m=1:2
    prior=old.reconstructedIncm{reference}(:,:,ix,m);
    for k=1:B
        coordinate=norm(lastAzel(:,:,k,m)-prior(:,:,k),'fro')/norm(prior(:,:,k),'fro');
        currentMetric=metric(lastAzel(:,:,k,m),k);priorMetric=metric(prior(:,:,k),k);
        db=abs(currentMetric.sinrDb-priorMetric.sinrDb);
        checks(end+1,:)={"coordinate_azel_surface",s.snrDb(k),measures(m),db,coordinate}; %#ok<AGROW>
        finer=extended(:,:,k,m);
        finerMetric=metric(finer,k);
        checks(end+1,:)={"surface_6401_12801",s.snrDb(k),measures(m), ...
            abs(finerMetric.sinrDb-priorMetric.sinrDb), ...
            norm(finer-prior(:,:,k),'fro')/norm(finer,'fro')}; %#ok<AGROW>
    end
    sectorIndex=find(string(labels(:,1))=="sector" & cell2mat(labels(:,3))==s.sectorVerificationHalfWidthDeg,1);
    coarse=records{sectorIndex};
    for k=1:B
        fineMetric=metric(sectorFine(:,:,k,m),k);coarseMetric=metric(coarse(:,:,k,m),k);
        checks(end+1,:)={"sector_3201_6401",s.snrDb(k),measures(m), ...
            abs(fineMetric.sinrDb-coarseMetric.sinrDb), ...
            norm(sectorFine(:,:,k,m)-coarse(:,:,k,m),'fro')/norm(sectorFine(:,:,k,m),'fro')}; %#ok<AGROW>
    end
end
checks=cell2table(checks,'VariableNames',{'comparison','inputSnrDb','measure','absSinrDeltaDb','relativeMatrixDelta'});
checks.passed=checks.absSinrDeltaDb<=s.coordinateSinrToleranceDb & ...
    checks.relativeMatrixDelta<=s.coordinateMatrixTolerance;
result=struct('cfg',c,'settings',s,'sourcePath',s.sourceJacobianResult, ...
    'inputCovariances',X,'truth',{truth},'records',{records},'labels',{labels}, ...
    'summary',summary,'powerTable',powerTable,'ablations',ablations,'contrast',contrast,'checks',checks, ...
    'singleSourceCovariances',single,'singleSourceReconstructions',singleRec, ...
    'sourceSnrDb',s.snrDb,'singleSourceQuadrature',singleQuadrature);
result.metadata=uav.experiment.metadata(c,"diagnosis_upa_integral_mechanism", ...
    struct('populationSnrPoints',B,'monteCarloTrials',0),toc(timer));
result.artifacts=struct('mat',fullfile(c.output.rawDir,[name '.mat']), ...
    'summaryCsv',fullfile(c.output.tableDir,[name '_summary.csv']), ...
    'powerCsv',fullfile(c.output.tableDir,[name '_power.csv']), ...
    'ablationCsv',fullfile(c.output.tableDir,[name '_ablation.csv']), ...
    'contrastCsv',fullfile(c.output.tableDir,[name '_contrast.csv']), ...
    'checksCsv',fullfile(c.output.tableDir,[name '_checks.csv']));
save(result.artifacts.mat,'result','-v7.3');
writetable(summary,result.artifacts.summaryCsv);writetable(powerTable,result.artifacts.powerCsv);
writetable(ablations,result.artifacts.ablationCsv);writetable(checks,result.artifacts.checksCsv);
writetable(contrast,result.artifacts.contrastCsv);
disp(checks);disp(powerTable);disp(ablations);disp(result.artifacts);
    function value=metric(R,k)
        [w,~]=uav.beamforming.solve("mvdr",R,c.derived.nominalSteeringReference,c);
        value=uav.metrics.evaluate(w,truth{k},c);
        value.noiseOutputPower=c.derived.noisePower*real(w'*w);
    end
    function appendAblation(label,k,measure,R)
        value=metric(R,k);
        ablationRows(end+1,:)={label,s.snrDb(k),measure,value.sinrDb, ...
            value.interferenceLeakage,value.noiseOutputPower};
    end
    function store(R,d,group,width)
        if nargin<4,width=c.spectrum.desiredExclusionHalfWidthDeg;end
        records{end+1}=R;labels(end+1,:)={group,d(1).gridSize,width};
        for inputIndex=1:size(X,3)
            if inputIndex<=B
                sourceIndex=inputIndex;label="full_population";
            elseif inputIndex==B+1
                sourceIndex=find(s.snrDb==10,1);label="true_incm_input";
            else
                sourceIndex=find(s.snrDb==10,1);label="white_input";
            end
            for measureIndex=1:numel(d)
                covariance=R(:,:,inputIndex,measureIndex);value=metric(covariance,sourceIndex);e=sort(real(eig(covariance)));
                shape=norm(covariance/trace(covariance)-truth{sourceIndex}.trueIncm/trace(truth{sourceIndex}.trueIncm),'fro');
                rows(end+1,:)={label,s.snrDb(sourceIndex),group,d(measureIndex).gridSize,width,d(measureIndex).measure, ...
                    value.sinrDb,value.interferenceLeakage,value.noiseOutputPower,e(end),mean(e(1:M-c.derived.J)),shape}; %#ok<AGROW>
            end
        end
        checkpoint=struct('cfg',c,'settings',s,'inputCovariances',X,'truth',{truth}, ...
            'records',{records},'labels',{labels},'summaryRows',{rows}, ...
            'generatedAt',string(datetime('now','TimeZone','Asia/Shanghai')));
        save(checkpointPath,'checkpoint','-v7.3');
        fprintf('原理诊断：%s，G=%d，排除半径=%g度，累计%.1f秒\n',group,d(1).gridSize,width,toc(timer));
    end
end
