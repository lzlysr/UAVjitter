function result=runRecComprehensiveCase(base,geo,condition)
%RUNRECCOMPREHENSIVECASE 独立条件块；所有方法、L和功率条件使用配对抽样。
g=base.experiments.recComprehensive; methods=g.methodNames; M=g.sensorCount;
cfg=base;
if condition.arrayType=="ULA", cfg.array.Nx=1; cfg.array.Ny=M; else, cfg.array.Nx=3; cfg.array.Ny=3; end
J=numel(condition.inrDb);
cfg.sources.desiredAzElDeg=geo.worldAzElDeg(1,:);
cfg.sources.interfererAzElDeg=geo.worldAzElDeg(2:J+1,:);
cfg.sources.interfererInrDb=condition.inrDb;
cfg.signal.snrDb=condition.snrDb;
cfg.attitude.stdDeg(:)=0;
cfg.timeline.L=max(g.snapshotCounts); % 固定全局最大长度，使跨条件的前缀抽样也一致。
cfg.spectrum.gridSize=condition.gridSize;
if condition.arrayType=="ULA", cfg.spectrum.gridSize=g.upaGridSize; end
cfg.spectrum.desiredExclusionHalfWidthDeg=condition.halfWidthDeg;
cfg.experiments.recComprehensive.ulaGridSize=condition.gridSize;
cfg=uav.config.validate(cfg);
gridCfg=cfg;
gridCfg.literature.guLeshem2012.upaDiagnosis.loadingRatio=gCandidate(base,'loadingRatio');
gridCfg.literature.guLeshem2012.upaDiagnosis.fitRegularization=gCandidate(base,'fitRegularization');
gridCfg.literature.guLeshem2012.upaDiagnosis.peakNmsRadiusBins=gCandidate(base,'peakNmsRadiusBins');
gridCfg.literature.guLeshem2012.upaDiagnosis.refineCandidates=false;
gridCfg.experiments.recUlaGeometry.candidateGridSize=condition.gridSize;
gridCfg.experiments.recUlaGeometry.loadingRatio=gCandidate(base,'loadingRatio');
gridCfg.experiments.recUlaGeometry.fitRegularization=gCandidate(base,'fitRegularization');
gridCfg.experiments.recUlaGeometry.peakNmsRadiusBins=gCandidate(base,'peakNmsRadiusBins');
historical=cfg;
historical.spectrum.gridSize=base.literature.guLeshem2012.upaDiagnosis.selectedGridSize;
historical=uav.config.validate(historical);
historical.experiments.recUlaGeometry.candidateGridSize=condition.gridSize;
historical.experiments.recUlaGeometry.loadingRatio=base.literature.guLeshem2012.upaDiagnosis.loadingRatio;
historical.experiments.recUlaGeometry.fitRegularization=base.literature.guLeshem2012.upaDiagnosis.fitRegularization;
historical.experiments.recUlaGeometry.peakNmsRadiusBins=base.literature.guLeshem2012.upaDiagnosis.peakNmsRadiusBins;
historical.literature.guLeshem2012.upaDiagnosis.refineCandidates=false;
if condition.arrayType=="UPA"
    rawGrid=uav.spectrum.upaIntegrationGrid(cfg,condition.gridSize,"solid_angle");
    scan=rawGrid.manifold(:,rawGrid.includedMask); weights=rawGrid.weights(rawGrid.includedMask);
else
    step=pi/condition.gridSize; angles=-pi/2+step/2+step*(0:condition.gridSize-1);
    reference=asin(cfg.derived.nominalDirectionCosines(2,1));
    keep=abs(angles-reference)>deg2rad(condition.halfWidthDeg);
    scan=uav.array.steeringULA(angles(keep),M,cfg.array.spacingOverLambda(2));
    weights=step*ones(1,sum(keep));
end
dimensions=[condition.trials,numel(condition.snapshotCounts),numel(methods)];
names=["sinrLinear","oracleLossDb","desiredGain","interferenceLeakage","noiseOutput", ...
    "incmNmse","shapeNmse","minimumEigenvalue","conditionNumber","weightLoading","fitKkt","runtimeSec"];
for name=names, result.(name)=nan(dimensions); end
result.status=zeros(dimensions,'int8'); % 1成功，0原式不适用，-1数值失败。
result.errors=cell(dimensions);
result.oracleSinrLinear=nan(condition.trials,1);
result.trueIncmNominalSinrLinear=nan(condition.trials,1);
result.desiredOutside=false(condition.trials,1); result.interfererInside=false(condition.trials,1);
result.attitudeRad=zeros(3,condition.trials);
nominal=partition(cfg.derived.nominalDirectionCosines,cfg,condition.arrayType);
% partition返回“在期望区外”；期望应在内，干扰应在外。
result.nominalPriorValid=~nominal(1) && all(nominal(2:end));
result.maxDesiredCoherence=0;
if J>0
    result.maxDesiredCoherence=max(abs(cfg.derived.nominalSteering(:,1)'*cfg.derived.nominalSteering(:,2:end)).^2/M^2);
end
result.minimumSourceGramEigenvalue=min(real(eig(cfg.derived.nominalSteering'*cfg.derived.nominalSteering/M)));
timer=tic;
for trial=1:condition.trials
    trialId=trial; if condition.population, trialId=0; end
    offset=g.snapshotSeedOffset+geo.id*g.geometryStride+trialId*g.trialStride;
    attitudeStream=uav.util.makeStream(cfg.random.baseSeed,offset+2);
    eta=deg2rad(condition.jitterStdDeg)*randn(attitudeStream,3,1);
    rotation=uav.geometry.rotationZYX(deg2rad(cfg.attitude.nominalDeg(:))+eta);
    u=rotation.'*cfg.derived.worldUnitVectors;
    [xi,directions]=uav.geometry.directionCoordinates(u);
    A=uav.array.steeringUPA(xi,cfg.array.Nx,cfg.array.Ny,cfg.array.spacingOverLambda);
    scenario=struct('steeringTruth',A,'arrayUnitVectors',u,'directionCosines',xi,'azElRad',directions);
    obs=uav.scenario.generatePeriod(cfg,scenario,1,uav.util.makeStream(cfg.random.baseSeed,offset+1));
    Ri=obs.truth.trueIncm; a=A(:,1); ps=cfg.derived.desiredPower;
    result.oracleSinrLinear(trial)=ps*real(a'*(Ri\a));
    nominalWeight=uav.beamforming.solve('mvdr',Ri,cfg.derived.nominalSteeringReference,cfg);
    result.trueIncmNominalSinrLinear(trial)=ps*abs(nominalWeight'*a)^2/real(nominalWeight'*Ri*nominalWeight);
    outside=partition(xi,cfg,condition.arrayType);
    result.desiredOutside(trial)=outside(1); result.interfererInside(trial)=any(~outside(2:end));
    result.attitudeRad(:,trial)=eta;
    for li=1:numel(condition.snapshotCounts)
        L=condition.snapshotCounts(li); small=obs;
        if isinf(L)
            small.Rx=Ri+ps*(a*a');
        else
            small.X=obs.X(:,1:L); small.pilot=obs.pilot(:,1:L); small.Rx=small.X*small.X'/L;
        end
        safe=uav.experiment.algorithmObservation(small);
        for mi=1:numel(methods)
            if methods(mi)=="raw_rec" && L<M, continue; end
            methodTimer=tic; kkt=NaN;
            try
                switch methods(mi)
                    case "raw_rec"
                        Rh=uav.incm.guLeshem2012Rec(safe.Rx,scan,weights);
                    case "old_static"
                        if condition.arrayType=="UPA"
                            Rh=uav.incm.reconstruct('static',safe,[],cfg);
                        else
                            Rh=uav.incm.ulaStatic(safe,reference,cfg,condition.halfWidthDeg);
                        end
                    case {"grid_covfit","historical_covfit"}
                        point=gridCfg; if methods(mi)=="historical_covfit", point=historical; end
                        if condition.arrayType=="UPA"
                            [Rh,d]=uav.incm.upaCovfit(safe,point);
                            B=d.candidateManifold; powers=d.fittedPowers; noise=d.noisePower;
                            regularization=point.literature.guLeshem2012.upaDiagnosis.fitRegularization;
                        else
                            [Rh,d]=uav.incm.ulaCovfit(safe,reference,point,condition.halfWidthDeg,false);
                            B=uav.array.steeringULA(d.nodesRad,M,cfg.array.spacingOverLambda(2));
                            powers=d.powers; noise=d.noisePower;
                            regularization=point.experiments.recUlaGeometry.fitRegularization;
                        end
                        assert(~d.usedTruth,'非Oracle方法不应读取真值。');
                        H=abs(B'*B).^2+regularization*eye(size(B,2));
                        b=real(diag(B'*(safe.Rx-noise*eye(M))*B)); gradient=H*powers-b;
                        scale=max([norm(H,inf)*norm(powers,inf),norm(b,inf),eps]);
                        kkt=max([max([0;-gradient])/scale; ...
                            norm(powers.*gradient,inf)/(max(norm(powers,inf),eps)*scale)]);
                    otherwise
                        error('uav:experiment:UnknownMethod','未知的比较方法。');
                end
                [w,wd]=uav.beamforming.solve('mvdr',Rh,cfg.derived.nominalSteeringReference,cfg);
                sinr=ps*abs(w'*a)^2/real(w'*Ri*w);
                assert(isfinite(sinr) && sinr>0 && sinr<=result.oracleSinrLinear(trial)*(1+1e-7));
                if condition.jitterStdDeg==0, assert(abs(w'*a-1)<1e-7); end
                result.sinrLinear(trial,li,mi)=sinr;
                result.oracleLossDb(trial,li,mi)=10*log10(result.oracleSinrLinear(trial)/sinr);
                result.desiredGain(trial,li,mi)=abs(w'*a)^2;
                result.interferenceLeakage(trial,li,mi)=real(w'*(Ri-cfg.derived.noisePower*eye(M))*w);
                result.noiseOutput(trial,li,mi)=cfg.derived.noisePower*real(w'*w);
                result.incmNmse(trial,li,mi)=norm(Rh-Ri,'fro')^2/norm(Ri,'fro')^2;
                align=real(sum(conj(Rh).*Ri,'all'))/norm(Rh,'fro')^2;
                result.shapeNmse(trial,li,mi)=norm(align*Rh-Ri,'fro')^2/norm(Ri,'fro')^2;
                result.minimumEigenvalue(trial,li,mi)=min(real(eig(Rh)));
                result.conditionNumber(trial,li,mi)=cond(Rh);
                result.weightLoading(trial,li,mi)=wd.loading;
                result.fitKkt(trial,li,mi)=kkt; result.status(trial,li,mi)=1;
            catch exception
                result.status(trial,li,mi)=-1;
                result.errors{trial,li,mi}=struct('identifier',exception.identifier,'message',exception.message);
            end
            result.runtimeSec(trial,li,mi)=toc(methodTimer);
        end
    end
end
result.condition=condition; result.geometry=geo; result.methodNames=methods;
result.elapsedSec=toc(timer);
result.snapshotSubstreams=g.snapshotSeedOffset+geo.id*g.geometryStride+ ...
    (1:condition.trials)*g.trialStride+1;
if condition.population, result.snapshotSubstreams(:)=g.snapshotSeedOffset+geo.id*g.geometryStride+1; end
end

function outside=partition(xi,cfg,arrayType)
%PARTITION 只用于诊断先验，不把真实源的分组信息传回算法。
if arrayType=="ULA"
    theta=asin(max(-1,min(1,xi(2,:))));
    reference=asin(cfg.derived.nominalDirectionCosines(2,1));
    outside=abs(theta-reference)>deg2rad(cfg.spectrum.desiredExclusionHalfWidthDeg);
else
    outside=uav.spectrum.admissibleDirections(xi.',cfg).';
end
end

function value=gCandidate(cfg,field)
%GCANDIDATE 明确复用此前冻结的候选配置，不按新几何调参。
value=cfg.experiments.recRootCause.candidateCovfit.(field);
end
