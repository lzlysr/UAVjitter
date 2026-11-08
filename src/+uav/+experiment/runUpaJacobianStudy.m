function result=runUpaJacobianStudy(cfg)
%RUNUPAJACOBIANSTUDY 同节点、同观测比较方向余弦面积与立体角积分。
%   只读取上轮已收敛的协方差，不重新生成信号、姿态或更换几何。
%   dOmega=du dv/sqrt(1-u^2-v^2)=dphi dz；采用后者避免边缘奇异。
%   原测度在相同节点上的权重为|z|dphi dz，作为独立坐标一致性控制。
timer=tic;s=cfg.experiments.upaJacobianStudy;
sourcePath=fullfile(cfg.output.rootDir,s.sourceGridResult);
loaded=load(sourcePath,'result');base=loaded.result;c=base.cfg;
if c.array.Nx<=1 || c.array.Ny<=1 || ~base.allChecksPassed
    error('uav:experiment:InvalidJacobianSource','必须读取已通过收敛检查的UPA诊断。');
end
c.experiments.upaJacobianStudy=s;
kernel=cfg.experiments.upaGridConvergence;
measures=["dxi_x_dxi_y","solid_angle"];
X=base.inputCovariances;N=size(X,3);
stamp=char(datetime('now','TimeZone','Asia/Shanghai','Format','yyyyMMdd_HHmmss'));
name=['diagnosis_upa_jacobian_' stamp];
checkpointPath=fullfile(c.output.rawDir,[name '_checkpoint.mat']);
reference=find(base.rules=="endpoint",1,'last');
cartesian=base.reconstructedIncm(:,:,:,reference);
cartSinr=base.sinrLinear(:,reference);
grids=s.gridSizes;rules=repmat("surface_midpoint",size(grids));
Rlist={};metrics={};quadrature={};runtimes=[];
cached=[];reused=[];
if strlength(s.reuseResult)>0
    previous=load(s.reuseResult,'result');cached=previous.result;
    if ~isequal(cached.inputCovariances,X) || ~isequal(cached.measures,measures) || ...
            ~isequal(cached.cfg.array,c.array) || ...
            ~isequal(cached.cfg.derived.nominalSteeringReference,c.derived.nominalSteeringReference) || ...
            ~isequal(cached.cfg.derived.nominalArrayUnitVectors,c.derived.nominalArrayUnitVectors) || ...
            cached.cfg.spectrum.desiredExclusionHalfWidthDeg~=c.spectrum.desiredExclusionHalfWidthDeg
        error('uav:experiment:InvalidJacobianCache','复用结果的观测、阵列或积分域不一致。');
    end
end
previousThreads=maxNumCompThreads(kernel.threadCount);
restore=onCleanup(@()maxNumCompThreads(previousThreads));
% 先运行预定网格序列与最细半格错位，再按预定上限决定是否继续加密。
for gridIndex=1:numel(grids)
    compute(gridIndex);
end
appendShift();
checks=checkConvergence();
if ~all(checks.passed) && max(s.extensionGridSizes)>max(grids)
    fprintf('Jacobian初轮未全部通过，增加预设的高密度复核，不放宽阈值。\n');
    for next=s.extensionGridSizes
        if next<=max(grids),continue;end
        grids(end+1)=next;rules(end+1)="surface_midpoint";compute(numel(grids)); %#ok<AGROW>
    end
    appendShift();checks=checkConvergence();
end
last=find(rules=="surface_midpoint",1,'last');
[coordinateDb,coordinateMatrix]=delta(Rlist{last}(:,:,:,1),cartesian,metrics{last}.sinrLinear(:,1),cartSinr);
coordinateCheck=struct('maxAbsSinrDb',max(coordinateDb),'maxRelativeMatrix',max(coordinateMatrix), ...
    'sinrToleranceDb',s.coordinateSinrToleranceDb,'matrixTolerance',s.coordinateMatrixRelativeTolerance, ...
    'passed',max(coordinateDb)<=s.coordinateSinrToleranceDb && max(coordinateMatrix)<=s.coordinateMatrixRelativeTolerance);
% 独立用原始求积核心核对两测度的权重，避免代数加速掩盖实现错误。
G=s.parityGridSize;az=-pi+pi/G+(0:G-1)*2*pi/G;
z=-1+1/(2*G)+(0:G-1)/G;[az,z]=meshgrid(az,z);
radius=sqrt(1-z(:).^2);nodes=[radius.*cos(az(:)),radius.*sin(az(:))];
keep=uav.spectrum.admissibleDirections(nodes,c);nodes=nodes(keep,:);z=z(keep);
A=uav.array.steeringUPA(nodes.',c.array.Nx,c.array.Ny,c.array.spacingOverLambda);
weights=[abs(z),ones(size(z))]*2*pi/G^2;
fast=uav.incm.integrateUpaBatch(X(:,:,1),c,G,"surface_midpoint",kernel,measures);
parity=zeros(1,numel(measures));
for parityMeasure=1:numel(measures)
    direct=uav.incm.guLeshem2012Rec(X(:,:,1),A,weights(:,parityMeasure));
    parity(parityMeasure)=norm(fast(:,:,1,parityMeasure)-direct,'fro')/norm(direct,'fro');
end
if any(parity>s.parityTolerance),error('uav:experiment:JacobianParity','与直接积分的独立核对失败。');end
summary=makeSummary();comparison=makeComparison(last);
result=struct('cfg',c,'settings',s,'sourceGridPath',sourcePath,'sourceGridMetadata',base.metadata, ...
    'sampleIndex',base.sampleIndex,'inputCovariances',X,'truth',{base.truth}, ...
    'measures',measures,'gridSizes',grids,'rules',rules,'quadrature',{quadrature}, ...
    'reconstructedIncm',{Rlist},'metrics',{metrics},'cartesianReference',cartesian, ...
    'cartesianReferenceSinrLinear',cartSinr,'oracleLinear',base.oracleLinear, ...
    'trueIncmFixedReferenceLinear',base.trueIncmFixedReferenceLinear, ...
    'summary',summary,'comparison',comparison,'checks',checks,'coordinateCheck',coordinateCheck, ...
    'allChecksPassed',all(checks.passed) && coordinateCheck.passed, ...
    'reusedGridResults',reused, ...
    'productionParityRelativeError',parity,'checkpointPath',checkpointPath);
result.metadata=uav.experiment.metadata(c,"diagnosis_upa_jacobian",base.metadata.actualCounts,toc(timer));
result.artifacts=struct('mat',fullfile(c.output.rawDir,[name '.mat']), ...
    'summaryCsv',fullfile(c.output.tableDir,[name '_summary.csv']), ...
    'comparisonCsv',fullfile(c.output.tableDir,[name '_comparison.csv']), ...
    'checksCsv',fullfile(c.output.tableDir,[name '_checks.csv']));
save(result.artifacts.mat,'result','-v7.3');writetable(summary,result.artifacts.summaryCsv);
writetable(comparison,result.artifacts.comparisonCsv);writetable(checks,result.artifacts.checksCsv);
fprintf('Jacobian诊断完成，收敛=%d，坐标一致性=%d，公式一致性最大误差=%.3g\n', ...
    all(checks.passed),coordinateCheck.passed,max(parity));
disp(comparison);disp(result.artifacts);
    function compute(g)
        gridTimer=tic;
        match=[];
        if ~isempty(cached),match=find(cached.gridSizes==grids(g) & cached.rules==rules(g),1);end
        reused(g)=~isempty(match);
        if isempty(match)
            [Rlist{g},quadrature{g}]=uav.incm.integrateUpaBatch(X,c,grids(g),rules(g),kernel,measures);
        else
            Rlist{g}=cached.reconstructedIncm{match};quadrature{g}=cached.quadrature{match};
        end
        record=struct('sinrLinear',zeros(N,2),'interferenceLeakage',zeros(N,2), ...
            'desiredGain',zeros(N,2),'shapeError',zeros(N,2),'numericalLoading',zeros(N,2));
        for n=1:N
            for m=1:numel(measures)
                R=Rlist{g}(:,:,n,m);truth=base.truth{n};
                [w,d]=uav.beamforming.solve("mvdr",R,c.derived.nominalSteeringReference,c);
                value=uav.metrics.evaluate(w,truth,c);
                if value.sinrLinear>base.oracleLinear(n)*(1+c.experiments.oracleTolerance)
                    error('uav:experiment:OracleBoundViolation','诊断结果超过真实最优上界。');
                end
                record.sinrLinear(n,m)=value.sinrLinear;
                record.interferenceLeakage(n,m)=value.interferenceLeakage;
                record.desiredGain(n,m)=value.desiredGain;
                record.numericalLoading(n,m)=d.loading;
                record.shapeError(n,m)=norm(R/real(trace(R))-truth.trueIncm/real(trace(truth.trueIncm)),'fro');
            end
        end
        metrics{g}=record;runtimes(g)=toc(gridTimer);
        if reused(g)
            runtimes(g)=cached.summary.gridRuntimeSec(find(cached.summary.gridSize==grids(g) & cached.summary.rule==rules(g),1));
        end
        example=find(base.sampleIndex.kind=="nominal_population" & base.sampleIndex.snr==10,1);
        fprintf('Jacobian：%s G=%d，节点%d，耗时%.2f秒；零抖动总体SNR10：du dv=%.6f，dOmega=%.6f dB\n', ...
            rules(g),grids(g),quadrature{g}(1).includedNodeCount,runtimes(g), ...
            10*log10(record.sinrLinear(example,1)),10*log10(record.sinrLinear(example,2)));
        checkpoint=struct('status',"partial",'cfg',c,'settings',s,'sourcePath',sourcePath, ...
            'gridSizes',grids,'rules',rules,'completedGridCount',g,'reconstructedIncm',{Rlist}, ...
            'metrics',{metrics},'quadrature',{quadrature},'runtimes',runtimes, ...
            'generatedAt',string(datetime('now','TimeZone','Asia/Shanghai','Format','yyyy-MM-dd''T''HH:mm:ssXXX')));
        save(checkpointPath,'checkpoint','-v7.3');
    end
    function appendShift()
        lastGrid=find(rules=="surface_midpoint",1,'last');
        grids(end+1)=grids(lastGrid);rules(end+1)="surface_shifted";compute(numel(grids));
    end
    function checks=checkConvergence()
        regular=find(rules=="surface_midpoint");shift=find(rules=="surface_shifted",1,'last');
        pairs=[regular(end-2),regular(end-1);regular(end-1),regular(end);regular(end),shift];
        rows=cell(0,7);
        for p=1:size(pairs,1)
            a=pairs(p,1);b=pairs(p,2);
            for m=1:numel(measures)
                [db,rel]=delta(Rlist{a}(:,:,:,m),Rlist{b}(:,:,:,m),metrics{a}.sinrLinear(:,m),metrics{b}.sinrLinear(:,m));
                rows(end+1,:)={measures(m),grids(a),grids(b),rules(b),max(db),max(rel), ...
                    max(db)<=s.sinrToleranceDb && max(rel)<=s.matrixRelativeTolerance}; %#ok<AGROW>
            end
        end
        checks=cell2table(rows,'VariableNames',{'measure','fromGrid','toGrid','toRule', ...
            'maxAbsSinrDeltaDb','maxRelativeMatrixDelta','passed'});
    end
    function summary=makeSummary()
        rows=cell(0,13);
        for caseName=["sample","population","nominal_population"]
            for inputSnr=unique(base.sampleIndex.snr).'
                ix=base.sampleIndex.kind==caseName & base.sampleIndex.snr==inputSnr;
                for g=1:numel(grids)
                    for m=1:numel(measures)
                        values=metrics{g}.sinrLinear(ix,m);
                        rows(end+1,:)={caseName,inputSnr,measures(m),grids(g),rules(g),sum(ix), ...
                            10*log10(median(values)),10*log10(quantile(values,.05)), ...
                            median(metrics{g}.interferenceLeakage(ix,m)),median(metrics{g}.desiredGain(ix,m)), ...
                            median(metrics{g}.shapeError(ix,m)),quadrature{g}(m).includedArea,runtimes(g)}; %#ok<AGROW>
                    end
                end
            end
        end
        summary=cell2table(rows,'VariableNames',{'caseName','inputSnrDb','measure','gridSize','rule', ...
            'sampleCount','medianSinrDb','p05SinrDb','medianInterferenceLeakage','medianDesiredGain', ...
            'medianTraceNormalizedShapeError','includedArea','gridRuntimeSec'});
    end
    function comparison=makeComparison(last)
        rows=cell(0,13);
        for caseName=["sample","population","nominal_population"]
            for inputSnr=unique(base.sampleIndex.snr).'
                ix=base.sampleIndex.kind==caseName & base.sampleIndex.snr==inputSnr;
                xy=metrics{last}.sinrLinear(ix,1);omega=metrics{last}.sinrLinear(ix,2);
                rows(end+1,:)={caseName,inputSnr,sum(ix),10*log10(median(cartSinr(ix))), ...
                    10*log10(median(xy)),10*log10(median(omega)), ...
                    10*log10(median(base.oracleLinear(ix))),10*log10(median(base.trueIncmFixedReferenceLinear(ix))), ...
                    median(10*log10(omega./xy)),10*log10(quantile(xy,.05)),10*log10(quantile(omega,.05)), ...
                    median(metrics{last}.shapeError(ix,1)),median(metrics{last}.shapeError(ix,2))}; %#ok<AGROW>
            end
        end
        comparison=cell2table(rows,'VariableNames',{'caseName','inputSnrDb','sampleCount', ...
            'cartesianMedianDb','sameNodeDxiMedianDb','solidAngleMedianDb','oracleMedianDb', ...
            'trueIncmFixedReferenceMedianDb','medianPairedGainDb','dxiP05Db','solidAngleP05Db', ...
            'dxiMedianShapeError','solidAngleMedianShapeError'});
    end
end
function [db,relative]=delta(A,B,sinrA,sinrB)
% 矩阵相对差以第二个输入为参照；逐样本比较，不让中位数掩盖不收敛样本。
db=abs(10*log10(sinrA./sinrB));relative=zeros(size(A,3),1);
for n=1:numel(relative),relative(n)=norm(A(:,:,n)-B(:,:,n),'fro')/norm(B(:,:,n),'fro');end
end
