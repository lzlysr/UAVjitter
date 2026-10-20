function results = runUlaGeometryStudy(profile)
%RUNULAGEOMETRYSTUDY 检查ULA中Rec与点源CovFit的几何敏感性。
%   原式只做一维dtheta积分；CovFit为明确标记的一维移植，不冒称原论文。
arguments
    profile (1,1) string = "smoke"
end
cfg = config(profile); g = cfg.experiments.recUlaGeometry;
cfg.array.Nx = 1; cfg.array.Ny = g.sensorCount;
cfg.array.spacingOverLambda(:) = g.spacingOverLambda;
cfg.attitude.stdDeg(:) = 0;
cfg.timeline.L = max(g.snapshotCounts);
cfg = uav.config.validate(cfg);
geometries = makeGeometries(cfg,g);
records = struct([]); candidateDiagnostics = cell(numel(geometries),numel(g.halfWidthsDeg),2);
recordCount = 0;
expectedCount = numel(geometries)*numel(g.halfWidthsDeg)*(numel(g.rawGridSizes)+2+ ...
    g.trialCount*numel(g.snapshotCounts)*3);
root = fileparts(fileparts(mfilename('fullpath')));
directory = fullfile(root,g.outputDirectory,profile);
if ~isfolder(directory), mkdir(directory); end
timer = tic;
for geometryIndex = 1:numel(geometries)
    geo = geometries(geometryIndex);
    local = cfg;
    % 仅借用既有快拍生成器；其合同不允许高度角恰好为0。
    % 用近水平夹具并反算方位角，保持y轴方向余弦严格等于sin(theta)。
    azimuths = asind(sind(geo.anglesDeg)/cosd(g.fixtureElevationDeg));
    local.sources.desiredAzElDeg = [azimuths(1),g.fixtureElevationDeg];
    local.sources.interfererAzElDeg = [azimuths(2:end).',repmat(g.fixtureElevationDeg,2,1)];
    local = uav.config.validate(local);
    streamOffset = g.seedOffset;
    if geo.family=="matched_world_projection"
        % 阵型间也共享源符号与阵元噪声；不只保证ULA方法间配对。
        streamOffset = cfg.experiments.recGeometry.seedOffset;
    end
    scenario = uav.scenario.generate(local,uav.util.makeStream(cfg.random.baseSeed,streamOffset));
    A = uav.array.steeringULA(deg2rad(geo.anglesDeg),g.sensorCount,g.spacingOverLambda);
    assert(norm(A-local.derived.nominalSteering,'fro')<1e-12,'ULA与退化UPA流形不一致。');
    Ri = local.derived.noisePower*eye(g.sensorCount)+A(:,2:end)*diag(local.derived.interfererPowers)*A(:,2:end)';
    Rx = Ri+local.derived.desiredPower*A(:,1)*A(:,1)';
    oracle = local.derived.desiredPower*real(A(:,1)'*(Ri\A(:,1)));
    for trial = 0:g.trialCount
        obs = uav.scenario.generatePeriod(local,scenario,1, ...
            uav.util.makeStream(cfg.random.baseSeed,streamOffset+trial+1));
        if trial==0, lengths=Inf; else, lengths=g.snapshotCounts; end
        for L = lengths
            small=obs;
            if isinf(L)
                small.Rx=Rx;
            else
                small.X=obs.X(:,1:L); small.Rx=small.X*small.X'/L;
            end
            safe = uav.experiment.algorithmObservation(small);
            for halfWidth = g.halfWidthsDeg
                for method = ["raw_rec","covfit_grid","covfit_refined"]
                    if isinf(L) && method=="raw_rec", sizes=g.rawGridSizes; else, sizes=g.rawGridSizes(1); end
                    for gridSize=sizes
                        Rh=[]; status="ok"; d=struct;
                        if method=="raw_rec"
                            if L<g.sensorCount
                                status="not-applicable-L<M";
                            else
                                step=pi/gridSize;
                                theta=-pi/2+step/2+step*(0:gridSize-1);
                                keep=abs(rad2deg(theta)-geo.anglesDeg(1))>halfWidth;
                                G=uav.array.steeringULA(theta(keep),g.sensorCount,g.spacingOverLambda);
                                Rh=uav.incm.guLeshem2012Rec(safe.Rx,G,step*ones(1,sum(keep)));
                            end
                        else
                            [Rh,d]=uav.incm.ulaCovfit(safe,deg2rad(geo.anglesDeg(1)), ...
                                local,halfWidth,method=="covfit_refined");
                            assert(~d.usedTruth,'CovFit不得使用真值。');
                            if trial==0
                                candidateDiagnostics{geometryIndex,find(g.halfWidthsDeg==halfWidth), ...
                                    1+(method=="covfit_refined")}=d;
                            end
                            gridSize=0;
                        end
                        value=NaN; nmse=NaN; leakage=NaN; loading=NaN;
                        if ~isempty(Rh)
                            [w,wd]=uav.beamforming.solve("mvdr",Rh,A(:,1),local);
                            value=local.derived.desiredPower*abs(w'*A(:,1))^2/real(w'*Ri*w);
                            nmse=norm(Rh-Ri,'fro')^2/norm(Ri,'fro')^2;
                            leakage=real(w'*(Ri-local.derived.noisePower*eye(g.sensorCount))*w);
                            loading=wd.loading;
                            assert(isfinite(value) && value>0 && value<=oracle*(1+1e-8));
                            assert(abs(w'*A(:,1)-1)<1e-8,'准确参考应保持期望响应。');
                        end
                        record=struct('geometryIndex',geometryIndex,'family',geo.family, ...
                            'parameter',geo.parameter,'desiredDeg',geo.anglesDeg(1), ...
                            'interferer1Deg',geo.anglesDeg(2),'interferer2Deg',geo.anglesDeg(3), ...
                            'trial',trial,'L',L,'halfWidthDeg',halfWidth,'method',method, ...
                            'snapshotSubstream',streamOffset+trial+1, ...
                            'rawGridSize',gridSize,'status',status,'sinrLinear',value, ...
                            'sinrDb',10*log10(value),'oracleSinrDb',10*log10(oracle), ...
                            'oracleGapDb',10*log10(oracle/value),'incmNmse',nmse, ...
                            'interferenceLeakage',leakage,'weightLoading',loading, ...
                            'interfererInsideSector',any(abs(geo.anglesDeg(2:end)-geo.anglesDeg(1))<=halfWidth), ...
                            'interfererCoherence',abs(A(:,2)'*A(:,3))^2/g.sensorCount^2);
                        if isempty(records), records=repmat(record,expectedCount,1); end
                        recordCount=recordCount+1; records(recordCount)=record;
                    end
                end
            end
        end
    end
    fprintf('ULA_GEOMETRY %d/%d %s parameter=%g elapsed=%.1fs\n', ...
        geometryIndex,numel(geometries),geo.family,geo.parameter,toc(timer));
    checkpoint=struct('cfg',local,'records',struct2table(records(1:recordCount)),'geometries',geometries);
    save(fullfile(directory,'checkpoint.mat'),'checkpoint','-v7.3');
end
assert(recordCount==expectedCount,'记录数量与预设工作量不一致。');
results=struct('cfg',cfg,'geometries',geometries,'records',struct2table(records), ...
    'candidateDiagnostics',{candidateDiagnostics}, ...
    'metadata',uav.experiment.metadata(cfg,"ula-geometry",struct('trials',g.trialCount),toc(timer)), ...
    'experimentSource',fileread(mfilename('fullpath')+".m"), ...
    'algorithmSource',fileread(fullfile(root,'src','+uav','+incm','ulaCovfit.m')), ...
    'configurationSource',fileread(fullfile(root,'config.m')));
save(fullfile(directory,'ula_geometry.mat'),'results','-v7.3');
writetable(results.records,fullfile(directory,'ula_records.csv'));
summary=groupsummary(results.records,{'geometryIndex','family','parameter','L', ...
    'halfWidthDeg','method','rawGridSize'},'median','sinrLinear');
summary.medianSinrDb=10*log10(summary.median_sinrLinear);
writetable(summary,fullfile(directory,'ula_summary.csv'));
fprintf('ULA_OUTPUT %s\n',directory);
end

function geometries=makeGeometries(cfg,g)
%MAKEGEOMETRIES 预设一维扫描，并将上一轮22组世界来向投影到y轴ULA。
geometries=struct('family',{},'parameter',{},'anglesDeg',{});
for angle=g.desiredAngleSweepDeg, add("desired_angle",angle,[angle,g.baseInterferersDeg]); end
for angle=g.commonRotationDeg, add("common_angle_shift",angle,[g.baseDesiredDeg,g.baseInterferersDeg]+angle); end
for angle=g.interfererOffsetsDeg, add("interferer_offset",angle,[g.offsetDesiredDeg,angle,g.offsetOtherInterfererDeg]); end
for separation=g.interfererSeparationsDeg
    add("interferer_pair_separation",separation,[g.baseDesiredDeg,g.pairCenterDeg+[-separation,separation]/2]);
end
root=fileparts(fileparts(mfilename('fullpath')));
previous=load(fullfile(root,cfg.experiments.recGeometry.outputDirectory,'preliminary','geometry_study.mat'),'results');
for index=1:numel(previous.results.geometries)
    point=previous.results.geometries(index);
    directions=[point.desired;point.interferers];
    projected=asind(cosd(directions(:,2)).*sind(directions(:,1))).';
    add("matched_world_projection",index,projected);
end
    function add(family,parameter,angles)
        geometries(end+1)=struct('family',family,'parameter',parameter,'anglesDeg',angles);
    end
end
