function results = runRecGeometryStudy(profile)
%RUNRECGEOMETRYSTUDY 在固定UPA和功率下检查源几何对Rec及CovFit的影响。
%   零抖动、准确名义参考；精确总协方差是隔离采样误差的诊断输入。
%   有限快拍在几何之间使用相同源符号及噪声随机流，各L取同一观测前缀。
%   严格Rec在L<M时不适用，不使用加载或伪逆补齐。所有失败点保留。
arguments
    profile (1,1) string = "smoke"
end
cfg = config(profile);
g = cfg.experiments.recGeometry;
study = cfg.experiments.recRootCause;
cfg.attitude.stdDeg(:) = 0;
cfg.timeline.L = max(g.snapshotCounts);
cfg.spectrum.gridSize = cfg.literature.guLeshem2012.upaDiagnosis.selectedGridSize;
cfg = uav.config.validate(cfg);
geometries = buildGeometries(cfg, g, study);
methodNames = ["oracle", "raw_rec_25", "raw_rec_7", "old_static", ...
    "historical_covfit", "refined_covfit_25", "refined_covfit_7"];
records = struct([]);
timer = tic;
root = fileparts(fileparts(mfilename('fullpath')));
directory = fullfile(root, g.outputDirectory, profile);
if ~isfolder(directory), mkdir(directory); end
for geometryIndex = 1:numel(geometries)
    point = geometries(geometryIndex);
    local = cfg;
    local.sources.desiredAzElDeg = point.desired;
    local.sources.interfererAzElDeg = point.interferers;
    local = uav.config.validate(local);
    scenario = uav.scenario.generate(local, uav.util.makeStream(cfg.random.baseSeed, g.seedOffset));
    A = local.derived.nominalSteering;
    Ri = local.derived.noisePower*eye(local.derived.M) + ...
        A(:,2:end)*diag(local.derived.interfererPowers)*A(:,2:end)';
    Rx = Ri + local.derived.desiredPower*A(:,1)*A(:,1)';
    coherence = abs(A'*A).^2/local.derived.M^2;
    descriptors = struct('maxDesiredInterfererCoherence', max(coherence(1,2:end)), ...
        'interfererCoherence', coherence(2,3), ...
        'minDesiredInterfererAngleDeg', min(acosd(max(-1,min(1, ...
        local.derived.worldUnitVectors(:,1)'*local.derived.worldUnitVectors(:,2:end))))));
    grids = cell(numel(g.rawGridSizes),2);
    for gridIndex = 1:numel(g.rawGridSizes)
        for sectorIndex = 1:2
            gridCfg = local;
            if sectorIndex == 2
                gridCfg.spectrum.desiredExclusionHalfWidthDeg = study.narrowSectorHalfWidthDeg;
            end
            grids{gridIndex,sectorIndex} = uav.spectrum.upaIntegrationGrid( ...
                gridCfg, g.rawGridSizes(gridIndex), "solid_angle");
        end
    end
    % trial=0、L=Inf明确表示总体诊断，并非有限快拍算法获得了真值。
    for trial = 0:g.trialCount
        obs = uav.scenario.generatePeriod(local, scenario, 1, ...
            uav.util.makeStream(cfg.random.baseSeed, g.seedOffset+trial+1));
        if trial == 0, lengths = Inf; else, lengths = g.snapshotCounts; end
        for L = lengths
            small = obs;
            pointCfg = local;
            if isinf(L)
                small.Rx = Rx;
            else
                small.X = obs.X(:,1:L);
                small.pilot = obs.pilot(:,1:L);
                small.Rx = small.X*small.X'/L;
                pointCfg.timeline.L = L;
            end
            safe = uav.experiment.algorithmObservation(small);
            for method = methodNames
                if startsWith(method,"raw_rec") && isinf(L)
                    gridIndices = 1:numel(g.rawGridSizes);
                else
                    gridIndices = 1;
                end
                for gridIndex = gridIndices
                    Rhat = [];
                    status = "ok";
                    methodCfg = pointCfg;
                    methodTimer = tic;
                    if method == "oracle"
                        Rhat = Ri;
                    elseif startsWith(method,"raw_rec")
                        sectorIndex = 1 + (method == "raw_rec_7");
                        grid = grids{gridIndex,sectorIndex};
                        if ~isinf(L) && L < local.derived.M
                            status = "not-applicable-L<M";
                        else
                            mask = grid.includedMask;
                            Rhat = uav.incm.guLeshem2012Rec(safe.Rx, ...
                                grid.manifold(:,mask),grid.weights(mask));
                        end
                    elseif method == "old_static"
                        Rhat = uav.incm.reconstruct("static",safe,[],methodCfg);
                    else
                        if startsWith(method,"refined")
                            methodCfg.spectrum.gridSize = study.candidateCovfit.gridSize;
                            d = methodCfg.literature.guLeshem2012.upaDiagnosis;
                            d.loadingRatio = study.candidateCovfit.loadingRatio;
                            d.peakNmsRadiusBins = study.candidateCovfit.peakNmsRadiusBins;
                            d.fitRegularization = study.candidateCovfit.fitRegularization;
                            d.refineCandidates = true;
                            methodCfg.literature.guLeshem2012.upaDiagnosis = d;
                            if method == "refined_covfit_7"
                                methodCfg.spectrum.desiredExclusionHalfWidthDeg = study.narrowSectorHalfWidthDeg;
                            end
                            methodCfg = uav.config.validate(methodCfg);
                        end
                        [Rhat,diagnostics] = uav.incm.upaCovfit(safe,methodCfg);
                        assert(~diagnostics.usedTruth,'非Oracle算法不得读取真值。');
                    end
                    sinr = NaN; leak = NaN; noise = NaN; nmse = NaN; loading = NaN;
                    if ~isempty(Rhat)
                        [w,wd] = uav.beamforming.solve("mvdr",Rhat, ...
                            safe.nominalSteeringReference,methodCfg);
                        leak = real(w'*(Ri-local.derived.noisePower*eye(local.derived.M))*w);
                        noise = local.derived.noisePower*real(w'*w);
                        sinr = local.derived.desiredPower*abs(w'*A(:,1))^2/real(w'*Ri*w);
                        nmse = norm(Rhat-Ri,'fro')^2/norm(Ri,'fro')^2;
                        loading = wd.loading;
                        assert(isfinite(sinr) && sinr>0,'出现非有限SINR。');
                        oracleW = Ri\A(:,1);
                        oracleSinr = local.derived.desiredPower*real(A(:,1)'*oracleW);
                        assert(sinr<=oracleSinr*(1+1e-8),'结果超过真实Oracle上界。');
                        assert(abs(w'*A(:,1)-1)<1e-8,'零失配MVDR未保持期望响应。');
                    end
                    % 零表示该方法不使用这类网格，避免把历史101谱网格误标为201。
                    rawGridSize = 0; caponGridSize = 0;
                    if startsWith(method,"raw_rec")
                        rawGridSize = g.rawGridSizes(gridIndex);
                    elseif method ~= "oracle"
                        caponGridSize = methodCfg.spectrum.gridSize;
                    end
                    record = struct('geometryIndex',geometryIndex,'family',point.family, ...
                        'parameter',point.parameter,'trial',trial,'L',L,'method',method, ...
                        'rawGridSize',rawGridSize,'caponGridSize',caponGridSize,'status',status, ...
                        'sinrLinear',sinr,'sinrDb',10*log10(sinr), ...
                        'interferenceLeakage',leak,'noiseOutput',noise, ...
                        'incmNmse',nmse,'weightLoading',loading,'runtimeSec',toc(methodTimer), ...
                        'maxDesiredInterfererCoherence',descriptors.maxDesiredInterfererCoherence, ...
                        'interfererCoherence',descriptors.interfererCoherence, ...
                        'minDesiredInterfererAngleDeg',descriptors.minDesiredInterfererAngleDeg);
                    records = [records;record]; %#ok<AGROW>
                end
            end
        end
    end
    fprintf('GEOMETRY %d/%d family=%s parameter=%g elapsed=%.1fs\n', ...
        geometryIndex,numel(geometries),point.family,point.parameter,toc(timer));
    % 每个几何点都落盘，长实验中断也保留已经完成的数据。
    checkpoint = struct('cfg',cfg,'geometries',geometries,'records',struct2table(records));
    save(fullfile(directory,'checkpoint.mat'),'checkpoint','-v7.3');
end
results = struct('cfg',cfg,'geometries',geometries,'records',struct2table(records), ...
    'methodNames',methodNames,'metadata',uav.experiment.metadata(cfg,"geometry-study", ...
    struct('trials',g.trialCount,'geometryCount',numel(geometries)),toc(timer)), ...
    'experimentSource',fileread(mfilename('fullpath')+".m"), ...
    'configurationSource',fileread(fullfile(root,'config.m')));
save(fullfile(directory,'geometry_study.mat'),'results','-v7.3');
writetable(results.records,fullfile(directory,'geometry_records.csv'));
% 只在独立试次间取线性SINR中位数；不报告小样本显著性。
summary = groupsummary(results.records, ...
    {'geometryIndex','family','parameter','L','method','rawGridSize'},'median','sinrLinear');
summary.medianSinrDb = 10*log10(summary.median_sinrLinear);
writetable(summary,fullfile(directory,'geometry_summary.csv'));
fprintf('GEOMETRY_OUTPUT %s\n',directory);
end

function geometries = buildGeometries(cfg,g,study)
%BUILDGEOMETRIES 四种预设扫描，不按结果筛选；端点及重复基准均保留。
geometries = struct('family',{},'parameter',{},'desired',{},'interferers',{});
for t = g.interpolationFractions
    add("interpolation",t,(1-t)*cfg.sources.desiredAzElDeg+t*study.orthogonalDesiredAzElDeg, ...
        (1-t)*cfg.sources.interfererAzElDeg+t*study.orthogonalInterfererAzElDeg);
end
for angle = g.commonAzimuthRotationsDeg
    desired = cfg.sources.desiredAzElDeg; interferers = cfg.sources.interfererAzElDeg;
    desired(1) = desired(1)+angle; interferers(:,1) = interferers(:,1)+angle;
    add("common_azimuth_rotation",angle,desired,interferers);
end
for angle = g.firstInterfererAzimuthOffsetsDeg
    interferers = cfg.sources.interfererAzElDeg;
    interferers(1,1) = cfg.sources.desiredAzElDeg(1)+angle;
    add("first_interferer_azimuth",angle,cfg.sources.desiredAzElDeg,interferers);
end
for elevation = g.desiredElevationsDeg
    desired = cfg.sources.desiredAzElDeg; desired(2) = elevation;
    add("desired_elevation",elevation,desired,cfg.sources.interfererAzElDeg);
end
    function add(family,parameter,desired,interferers)
        geometries(end+1) = struct('family',family,'parameter',parameter, ...
            'desired',desired,'interferers',interferers);
    end
end
