function result=runUpaPopulationConvergence(cfg)
%RUNUPAPOPULATIONCONVERGENCE 同几何、零抖动、无采样误差的UPA积分收敛检查。
%   正式incm_integral_origin保持原式。此处只改变求积网格和测度；
%   不加载Capon、不修复重构矩阵、不拟合功率、不估计或修正期望方向。
%   方位/负z坐标与方位/高度角坐标独立求积，避免在地平线裁剪Jacobian。

timer=tic;s=cfg.experiments.upaPopulationConvergence;
c=cfg;c.array.Nx=s.arrayShape(1);c.array.Ny=s.arrayShape(2);
c.attitude.stdDeg=zeros(size(c.attitude.stdDeg));c=uav.config.validate(c);
if c.derived.arrayType~="UPA" || numel(s.gridSizes)<3 || ...
        any(diff(s.gridSizes)<=0) || any(s.gridSizes<3) || ...
        any(mod(s.gridSizes,2)~=1) || any(s.extensionGridSizes<=s.gridSizes(end)) || ...
        any(diff(s.extensionGridSizes)<=0)
    error('uav:experiment:InvalidPopulationStudy','需要二维UPA和至少三级递增奇数网格，扩展网格必须更密。');
end
% 用与场景模块相同的名义流形直接形成理论协方差；不生成随机快拍。
A=c.derived.nominalSteering;M=c.derived.M;B=numel(s.snrDb);
Ai=A(:,2:end);Rin=c.derived.noisePower*eye(M)+ ...
    Ai*diag(c.derived.interfererPowers)*Ai';
inputs=complex(zeros(M,M,B+2));truth=cell(1,B);oracle=zeros(1,B);
for b=1:B
    ps=c.derived.noisePower*10^(s.snrDb(b)/10);
    inputs(:,:,b)=Rin+ps*(A(:,1)*A(:,1)');
    truth{b}=struct('trueASVs',A,'trueIncm',Rin,'desiredPower',ps, ...
        'interfererPowers',c.derived.interfererPowers);
    w=uav.beamforming.solve("mvdr",Rin,A(:,1),c);
    v=uav.metrics.evaluate(w,truth{b},c);oracle(b)=v.sinrDb;
end
% 两个额外输入只作原理控制，不能标为实际算法的性能曲线。
inputs(:,:,B+1)=Rin;inputs(:,:,B+2)=c.derived.noisePower*eye(M);
root=fullfile(c.output.rootDir,s.outputSubdirectory);
if ~isfolder(root),mkdir(root);end
stamp=char(datetime('now','TimeZone','Asia/Shanghai','Format','yyyyMMdd_HHmmss_SSS'));
[~,suffix]=fileparts(tempname(root));
runDir=fullfile(root,sprintf('upa_%dx%d_population_%s_%s',c.array.Nx,c.array.Ny,stamp,suffix));
if isfolder(runDir),error('uav:experiment:DiagnosticCollision','独立结果目录已存在，拒绝覆盖。');end
mkdir(runDir);logPath=fullfile(runDir,'progress.log');
previousThreads=maxNumCompThreads(s.threadCount);
restoreThreads=onCleanup(@()maxNumCompThreads(previousThreads)); %#ok<NASGU>
records=struct('rule',{},'gridSize',{},'measures',{},'R',{},'quadrature',{}, ...
    'sinrDb',{},'leakage',{},'noiseOutputPower',{},'loading',{},'runtimeSec',{});
% 追加更密网格时可复用此前完成的求积，不重复消耗计算时间。
% 续算仍新建输出目录，且逐元素核对理论输入和排除区，拒绝混用场景。
if strlength(s.resumeCheckpoint)>0
    prior=load(s.resumeCheckpoint,'records');
    old=load(fullfile(fileparts(s.resumeCheckpoint),'inputs.mat'),'inputs','c');
    if ~isequal(inputs,old.inputs) || ...
            ~isequal(c.derived.nominalSteering,old.c.derived.nominalSteering) || ...
            c.spectrum.desiredExclusionHalfWidthDeg~=old.c.spectrum.desiredExclusionHalfWidthDeg
        error('uav:experiment:ResumeMismatch','续算输入或排除区与原检查点不一致。');
    end
    records=prior.records;
    for k=1:numel(records)
        if ~isfield(records(k).quadrature,'nativeAcceleration')
            [records(k).quadrature.nativeAcceleration]=deal(false);
        end
    end
end
say('开始：%dx%d，零抖动，理论协方差，%d个SNR点；输出%s', ...
    c.array.Nx,c.array.Ny,B,runDir);
save(fullfile(runDir,'inputs.mat'),'c','s','inputs','truth','oracle','-v7.3');
% 先与未经加速的原始生产公式逐元素核对；高功率输入最容易暴露消去误差。
parity=checkParity();
if any(parity.relativeMatrixError>s.parityTolerance)
    error('uav:experiment:PopulationParityFailure','分块加速未通过原始Rec公式一致性核对。');
end
writetable(parity,fullfile(runDir,'parity.csv'));
for g=s.gridSizes,compute(g,"endpoint","dxi_x_dxi_y");end
% 两种测度使用完全相同的surface节点；Jacobian通过|z|权重的有无体现。
for g=s.gridSizes
    compute(g,"surface_midpoint",["dxi_x_dxi_y","solid_angle"]);
end
% 未达阈值时追加精度；已经收敛的序列不必重复计算更密网格。
% 若坐标/错位检查未通过，也继续加密，而不是只检查一次便宣布完成。
for round=0:numel(s.extensionGridSizes)
    primary=convergenceChecks(false);
    if all(primary.passed) || round==numel(s.extensionGridSizes)
        surfaceIndex=find([records.rule]=="surface_midpoint",1,'last');
        finalGrid=records(surfaceIndex).gridSize;
        compute(finalGrid,"surface_shifted",["dxi_x_dxi_y","solid_angle"]);
        % 独立高度角求积使用两级；同样需通过矩阵及SINR阈值。
        compute((finalGrid-1)/2+1,"azel_midpoint",["dxi_x_dxi_y","solid_angle"]);
        compute(finalGrid,"azel_midpoint",["dxi_x_dxi_y","solid_angle"]);
        extra=convergenceChecks(true);coordinates=coordinateComparison();
        if all(extra.passed) && all(coordinates.passed),break;end
    end
    if round<numel(s.extensionGridSizes)
        g=s.extensionGridSizes(round+1);
        for rule=["endpoint","surface_midpoint"]
            needsPrimary=any(~primary.passed(primary.rule==rule));
            lastGrid=max([records([records.rule]==rule).gridSize]);
            if g>lastGrid && (needsPrimary || all(primary.passed))
                if rule=="endpoint",measures="dxi_x_dxi_y";
                else,measures=["dxi_x_dxi_y","solid_angle"];end
                compute(g,rule,measures);
            end
        end
    end
end
checks=convergenceChecks(true);coordinateChecks=coordinateComparison();
summary=makeSummary();finalComparison=makeFinalComparison();controls=makeControls();
allPassed=all(checks.passed) && all(coordinateChecks.passed);
result=struct('cfg',c,'settings',s,'inputCovariances',inputs,'truth',{truth}, ...
    'oracleSinrDb',oracle,'records',records,'parity',parity,'summary',summary, ...
    'checks',checks,'coordinateChecks',coordinateChecks,'comparison',finalComparison, ...
    'controls',controls,'allChecksPassed',allPassed,'outputDirectory',string(runDir));
result.metadata=uav.experiment.metadata(c,"upa_population_integral_convergence", ...
    struct('snrPoints',B,'monteCarloTrials',0,'jitterRmsDeg',0,'arrayElements',M),toc(timer));
result.metadata.cumulativeGridRuntimeSec=sum([records.runtimeSec]);
result.metadata.resumeCheckpoint=s.resumeCheckpoint;
save(fullfile(runDir,'study.mat'),'result','-v7.3');
writetable(summary,fullfile(runDir,'summary.csv'));
writetable(checks,fullfile(runDir,'convergence_checks.csv'));
writetable(coordinateChecks,fullfile(runDir,'coordinate_checks.csv'));
writetable(finalComparison,fullfile(runDir,'final_comparison.csv'));
writetable(controls,fullfile(runDir,'truth_controls.csv'));
drawFigures();
say('完成：数值检查全部通过=%d，耗时%.1f秒。',allPassed,toc(timer));
disp(finalComparison);disp(checks);disp(coordinateChecks);

    function compute(g,rule,measures)
        if ~isempty(records) && any([records.rule]==rule & [records.gridSize]==g),return;end
        stepTimer=tic;say('开始 %s，G=%d。',rule,g);
        [R,d]=uav.incm.integrateUpaBatch(inputs,c,g,rule,s,measures);
        rec=struct('rule',rule,'gridSize',g,'measures',measures,'R',R, ...
            'quadrature',d,'sinrDb',zeros(B,numel(measures)), ...
            'leakage',zeros(B,numel(measures)),'noiseOutputPower',zeros(B,numel(measures)), ...
            'loading',zeros(B,numel(measures)),'runtimeSec',toc(stepTimer));
        for b=1:B
            for m=1:numel(measures)
                [w,wd]=uav.beamforming.solve("mvdr",R(:,:,b,m),A(:,1),c);
                value=uav.metrics.evaluate(w,truth{b},c);
                rec.sinrDb(b,m)=value.sinrDb;rec.leakage(b,m)=value.interferenceLeakage;
                rec.noiseOutputPower(b,m)=c.derived.noisePower*real(w'*w);
                rec.loading(b,m)=wd.loading;
                if value.sinrDb>oracle(b)+c.experiments.oracleTolerance
                    error('uav:experiment:OracleBoundViolation','理论协方差诊断超过Oracle上界。');
                end
            end
        end
        records(end+1)=rec;
        % 只更新本次新建目录中的检查点，不触碰正式实验MAT或图形。
        save(fullfile(runDir,'checkpoint.mat'),'records','s','-v7.3');
        writetable(makeSummary(),fullfile(runDir,'summary.csv'));
        center=find(s.snrDb==10,1);if isempty(center),center=1;end
        say('完成 %s，G=%d，%.2f秒，SNR=%g：%s dB。', ...
            rule,g,rec.runtimeSec,s.snrDb(center),mat2str(rec.sinrDb(center,:),8));
    end
    function out=checkParity()
        rows={};G=s.parityGridSize;
        % 检查最高SNR输入及白噪声输入；同时覆盖端点与两种球面积分坐标。
        for rule=["endpoint","surface_midpoint","azel_midpoint"]
            if rule=="endpoint",measures="dxi_x_dxi_y";
            else,measures=["dxi_x_dxi_y","solid_angle"];end
            if rule=="endpoint"
                axis=linspace(-1,1,G);[x,y]=meshgrid(axis);nodes=[x(:),y(:)];
                weights=repmat((2/(G-1))^2,G^2,1);
            else
                az=-pi+pi/G+(0:G-1)*2*pi/G;
                if rule=="surface_midpoint"
                    z=-1+1/(2*G)+(0:G-1)/G;[phi,z]=meshgrid(az,z);
                    rr=sqrt(1-z(:).^2);nodes=[rr.*cos(phi(:)),rr.*sin(phi(:))];
                    weights=[abs(z(:)),ones(G^2,1)]*2*pi/G^2;
                else
                    el=-pi/2+pi/(4*G)+(0:G-1)*pi/(2*G);[phi,el]=meshgrid(az,el);
                    rr=cos(el(:));nodes=[rr.*cos(phi(:)),rr.*sin(phi(:))];
                    weights=[abs(sin(el(:))).*rr,rr]*pi^2/G^2;
                end
            end
            keep=uav.spectrum.admissibleDirections(nodes,c);nodes=nodes(keep,:);weights=weights(keep,:);
            manifold=uav.array.steeringUPA(nodes.',c.array.Nx,c.array.Ny,c.array.spacingOverLambda);
            for inputIndex=[B,B+2]
                fast=uav.incm.integrateUpaBatch(inputs(:,:,inputIndex),c,G,rule,s,measures);
                for m=1:numel(measures)
                    direct=uav.incm.guLeshem2012Rec(inputs(:,:,inputIndex),manifold,weights(:,m));
                    rel=norm(fast(:,:,1,m)-direct,'fro')/norm(direct,'fro');
                    rows(end+1,:)={rule,measures(m),inputIndex,rel}; %#ok<AGROW>
                end
            end
        end
        out=cell2table(rows,'VariableNames',{'rule','measure','inputIndex','relativeMatrixError'});
        say('原式一致性核对：最大相对误差%.3g。',max(out.relativeMatrixError));
    end
    function out=convergenceChecks(includeExtra)
        rows={};
        for rule=["endpoint","surface_midpoint"]
            ix=find([records.rule]==rule);
            if numel(ix)<3,continue;end
            for p=numel(ix)-2:numel(ix)-1
                append(ix(p),ix(p+1),"successive_refinement");
            end
        end
        if includeExtra
            last=find([records.rule]=="surface_midpoint",1,'last');
            shifted=find([records.rule]=="surface_shifted",1,'last');
            append(last,shifted,"half_azimuth_shift");
            fine=records(last).gridSize;coarse=(fine-1)/2+1;
            i=find([records.rule]=="azel_midpoint" & [records.gridSize]==coarse,1);
            j=find([records.rule]=="azel_midpoint" & [records.gridSize]==fine,1);
            append(i,j,"independent_azel_refinement");
        end
        out=cell2table(rows,'VariableNames',{'comparison','rule','measure','fromGrid','toGrid', ...
            'maxSinrDeltaDb','maxRelativeMatrixDelta','passed'});
        function append(i,j,label)
            a=records(i);b=records(j);
            for m=1:numel(a.measures)
                deltaDb=max(abs(a.sinrDb(:,m)-b.sinrDb(:,m)));
                deltaR=matrixDelta(a.R(:,:,:,m),b.R(:,:,:,m));
                rows(end+1,:)={label,b.rule,a.measures(m),a.gridSize,b.gridSize,deltaDb,deltaR, ...
                    deltaDb<=s.sinrToleranceDb && deltaR<=s.matrixRelativeTolerance}; %#ok<AGROW>
            end
        end
    end
    function out=coordinateComparison()
        e=find([records.rule]=="endpoint",1,'last');
        f=find([records.rule]=="surface_midpoint",1,'last');
        a=find([records.rule]=="azel_midpoint",1,'last');rows={};
        pairs=[e,f,1;f,a,1;f,a,2];
        for k=1:size(pairs,1)
            i=pairs(k,1);j=pairs(k,2);m=pairs(k,3);
            db=max(abs(records(i).sinrDb(:,m)-records(j).sinrDb(:,m)));
            deltaR=matrixDelta(records(i).R(:,:,:,m),records(j).R(:,:,:,m));
            rows(end+1,:)={records(i).rule,records(j).rule,records(j).measures(m),db,deltaR, ...
                db<=s.coordinateSinrToleranceDb && deltaR<=s.coordinateMatrixTolerance}; %#ok<AGROW>
        end
        % 球冠完全位于下半球时，补域面积和立体角有解析值，可独立检查权重。
        beta=deg2rad(c.spectrum.desiredExclusionHalfWidthDeg);
        zs=abs(c.derived.nominalArrayUnitVectors(3,1));
        if asin(zs)<=beta
            error('uav:experiment:CapCrossesHorizon','本解析域检查要求期望排除球冠完全位于下半球。');
        end
        exact=[pi-pi*sin(beta)^2*zs,2*pi*cos(beta)];
        for m=1:2
            rel=abs(records(f).quadrature(m).includedArea/exact(m)-1);
            rows(end+1,:)={"analytic_domain",records(f).rule,records(f).measures(m),NaN,rel, ...
                rel<=s.domainRelativeTolerance}; %#ok<AGROW>
        end
        out=cell2table(rows,'VariableNames',{'fromRule','toRule','measure','maxSinrDeltaDb', ...
            'maxRelativeMatrixOrDomainError','passed'});
    end
    function out=makeSummary()
        rows={};
        for k=1:numel(records)
            rec=records(k);
            for b=1:B
                for m=1:numel(rec.measures)
                    R=rec.R(:,:,b,m);shape=norm(R/real(trace(R))-Rin/real(trace(Rin)),'fro');
                    rows(end+1,:)={rec.rule,rec.gridSize,rec.measures(m),s.snrDb(b),rec.sinrDb(b,m), ...
                        oracle(b),oracle(b)-rec.sinrDb(b,m),rec.leakage(b,m),rec.noiseOutputPower(b,m), ...
                        rec.loading(b,m),shape,rec.quadrature(m).includedArea,rec.runtimeSec}; %#ok<AGROW>
                end
            end
        end
        out=cell2table(rows,'VariableNames',{'rule','gridSize','measure','inputSnrDb','sinrDb', ...
            'oracleSinrDb','gapToOracleDb','interferenceLeakage','noiseOutputPower','numericalLoading', ...
            'traceNormalizedShapeError','domainMass','gridRuntimeSec'});
    end
    function out=makeFinalComparison()
        e=find([records.rule]=="endpoint",1,'last');f=find([records.rule]=="surface_midpoint",1,'last');
        out=table(s.snrDb(:),oracle(:),records(e).sinrDb(:,1),records(f).sinrDb(:,1), ...
            records(f).sinrDb(:,2),records(f).sinrDb(:,2)-records(f).sinrDb(:,1), ...
            'VariableNames',{'inputSnrDb','oracleSinrDb','originalAreaSinrDb','sameNodesAreaSinrDb', ...
            'solidAngleSinrDb','jacobianGainDb'});
    end
    function out=makeControls()
        % 仅用于诊断的真值对照，绝不能作为无先验的实际方法参与排名。
        f=find([records.rule]=="surface_midpoint",1,'last');rec=records(f);
        rows=zeros(B,6);
        for b=1:B
            denominator=real(sum(conj(Ai).*(inputs(:,:,b)\Ai),1));
            pointPower=1./denominator;
            pointR=c.derived.noisePower*eye(M)+Ai*diag(pointPower)*Ai';
            w=uav.beamforming.solve("mvdr",pointR,A(:,1),c);
            v=uav.metrics.evaluate(w,truth{b},c);
            rows(b,1:4)=[s.snrDb(b),oracle(b),v.sinrDb, ...
                max(abs(pointPower./c.derived.interfererPowers(:)'-1))];
            for m=1:2
                % 输入已经移除期望源的理论INCM，仍按同一积分公式重构。
                w=uav.beamforming.solve("mvdr",rec.R(:,:,B+1,m),A(:,1),c);
                v=uav.metrics.evaluate(w,truth{b},c);rows(b,4+m)=v.sinrDb;
            end
        end
        out=array2table(rows,'VariableNames',{'inputSnrDb','oracleSinrDb', ...
            'trueDirectionPointCaponSinrDb','maxPointPowerRelativeError', ...
            'noDesiredInputAreaSinrDb','noDesiredInputSolidAngleSinrDb'});
    end
    function drawFigures()
        pc=c;pc.plot.visible=false;pc.plot.closeAfterExport=true;
        pc.plot.figurePosition=s.figurePosition;pc.plot.fontSize=s.figureFontSize;
        fig=uav.plot.newFigure(pc);hold on;
        center=find(s.snrDb==10,1);if isempty(center),center=1;end
        styles={'-s','--o',':^','-.d','-x'};series=0;
        for rule=["endpoint","surface_midpoint","azel_midpoint"]
            ix=find([records.rule]==rule);gs=[records(ix).gridSize];
            [gs,order]=sort(gs);ix=ix(order);
            for m=1:numel(records(ix(1)).measures)
                ys=arrayfun(@(i)records(i).sinrDb(center,m),ix);
                if rule=="endpoint",label="方向余弦坐标：面积积分";
                elseif rule=="surface_midpoint",label="方位—负z坐标：";
                else,label="方位—高度角坐标：";end
                if rule~="endpoint"
                    if records(ix(1)).measures(m)=="solid_angle",label=label+"立体角积分";
                    else,label=label+"面积积分";end
                end
                series=series+1;
                semilogx(gs,ys,styles{series},'DisplayName',label,'LineWidth',c.plot.lineWidth);
            end
        end
        yline(oracle(center),'k--','DisplayName','Oracle（理想上界）');
        xlabel('每轴求积点数G（对数坐标）');ylabel('输出SINR（dB）');grid on;
        set(gca,'FontSize',s.figureFontSize,'XScale','log');legend('Location','eastoutside','FontSize',s.figureFontSize,'Interpreter','none');
        status="尚未全部通过";if allPassed,status="全部数值检查通过";end
        title(sprintf('零抖动、理论协方差，SNR=%g dB；%s',s.snrDb(center),status));
        export(fig,'grid_convergence');
        fig=uav.plot.newFigure(pc);hold on;
        plot(s.snrDb,oracle,'k--','LineWidth',c.plot.lineWidth,'DisplayName','Oracle（理想上界）');
        labels={'原始方向余弦面积积分','同节点面积积分','立体角积分（Jacobian）'};
        styles={'--s',':o','-.^'};
        for k=1:3
            plot(finalComparison.inputSnrDb,finalComparison{:,k+2},styles{k}, ...
                'LineWidth',c.plot.lineWidth,'DisplayName',labels{k});
        end
        xlabel('输入SNR（dB）');ylabel('输出SINR（dB）');grid on;
        set(gca,'FontSize',s.figureFontSize);legend('Location','eastoutside','FontSize',s.figureFontSize);
        title(sprintf('最高网格结果；%s',status));export(fig,'sinr_comparison');
    end
    function export(fig,name)
        cleanup=onCleanup(@()close(fig)); %#ok<NASGU>
        savefig(fig,fullfile(runDir,[name '.fig']));
        exportgraphics(fig,fullfile(runDir,[name '.pdf']),'ContentType','vector');
        exportgraphics(fig,fullfile(runDir,[name '.png']),'Resolution',c.plot.resolutionDpi);
    end
    function say(format,varargin)
        message=sprintf(format,varargin{:});
        when=char(datetime('now','TimeZone','Asia/Shanghai','Format','yyyy-MM-dd HH:mm:ss'));
        fprintf('[%s] %s\n',when,message);
        fid=fopen(logPath,'a');if fid<0,error('uav:experiment:LogWrite','无法写入诊断日志。');end
        cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
        fprintf(fid,'[%s] %s\n',when,message);
    end
end

function delta=matrixDelta(a,b)
% 同测度比较未归一化矩阵，防止只看SINR漏掉整体数值积分未收敛。
delta=0;
for k=1:size(a,3)
    delta=max(delta,norm(a(:,:,k)-b(:,:,k),'fro')/norm(b(:,:,k),'fro'));
end
end
