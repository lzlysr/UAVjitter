function result = runCovfitDirectionOnlyControl(resultFile)
%RUNCOVFITDIRECTIONONLYCONTROL 固定拟合功率，仅细化干扰候选方向的诊断。
%   每种L选择零抖动下“细化前后SINR差最大”的已记录点以放大症状；
%   这不是无偏性能统计，不能用所选点替代完整配对验证。
arguments
    resultFile (1,1) string
end
saved = load(resultFile,'results');
r = saved.results;
cfg = r.cfg;
s = cfg.experiments.recRootCause;
point = find(s.validationJitterDeg==0);
beforeIndex = find(r.methodNames=="candidate_covfit");
afterIndex = find(r.methodNames=="refined_candidate_covfit");
rows = struct([]);
for snapshotIndex = 1:numel(s.validationSnapshotCounts)
    before = r.sinrLinear(:,:,point,snapshotIndex,beforeIndex);
    after = r.sinrLinear(:,:,point,snapshotIndex,afterIndex);
    [~,index] = max(log(after(:)./before(:)));
    [trial,k] = ind2sub(size(before),index);
    local = cfg; local.attitude.stdDeg(:) = 0;
    scene = uav.scenario.generate(local,uav.util.makeStream(cfg.random.baseSeed,s.validationSeedOffset+trial));
    obs = uav.scenario.generatePeriod(local,scene,k,uav.util.makeStream( ...
        cfg.random.baseSeed,s.validationSeedOffset+s.validationTrials+(trial-1)*s.validationPeriods+k));
    L = s.validationSnapshotCounts(snapshotIndex);
    obs.X = obs.X(:,1:L); obs.Rx = obs.X*obs.X'/L; obs.pilot = obs.pilot(:,1:L);
    candidate = r.candidateCfg;
    candidate.timeline.L = L;
    candidate.literature.guLeshem2012.upaDiagnosis.refineCandidates = false;
    [Rh,d] = uav.incm.upaCovfit(uav.experiment.algorithmObservation(obs),candidate);
    loading = s.candidateCovfit.loadingRatio*real(trace(obs.Rx))/cfg.derived.M;
    [nodes,refinement] = uav.spectrum.refineCaponNodes(d.interfererCandidateNodes, ...
        obs.Rx+loading*eye(cfg.derived.M),false,candidate);
    Ai = uav.array.steeringUPA(nodes.',cfg.array.Nx,cfg.array.Ny,cfg.array.spacingOverLambda);
    % 候选数、逐候选系数、噪声、期望标签和权值参考全部保持不变。
    powers = d.interfererFittedPowers;
    directionOnly = (Ai.*powers.')*Ai'+d.noisePower*eye(cfg.derived.M);
    w = uav.beamforming.solve("mvdr",Rh,cfg.derived.nominalSteeringReference,cfg);
    m = uav.metrics.evaluate(w,obs.truth,cfg);
    wd = uav.beamforming.solve("mvdr",directionOnly,cfg.derived.nominalSteeringReference,cfg);
    md = uav.metrics.evaluate(wd,obs.truth,cfg);
    % 同时检查细化前NNLS的KKT，防止把优化器失效误称为方向误差。
    A = d.candidateManifold; q = d.fittedPowers;
    residual = (A.*q.')*A'+d.noisePower*eye(cfg.derived.M)-obs.Rx;
    gradient = 2*real(sum(conj(A).*(residual*A),1)).';
    target = obs.Rx-d.noisePower*eye(cfg.derived.M);
    scale = max(norm(2*real(sum(conj(A).*(target*A),1)),inf),eps);
    active = q>sqrt(eps)*max(norm(q,inf),1);
    kkt = max([max(-gradient,0);abs(gradient(active))])/scale;
    row = struct('L',L,'trial',trial,'period',k,'beforeSinrDb',m.sinrDb, ...
        'directionOnlySinrDb',md.sinrDb,'fullyRefinedSinrDb',10*log10(after(trial,k)), ...
        'beforeFitRelativeKktViolation',kkt,'interfererCandidateCount',numel(powers), ...
        'refinementAcceptedCount',sum(refinement.accepted), ...
        'beforeReplayRelativeError',abs(m.sinrLinear-before(trial,k))/before(trial,k));
    rows = [rows;row]; %#ok<AGROW>
end
result = struct('sourceValidation',resultFile,'cfg',cfg,'records',struct2table(rows), ...
    'selectionRule',"每种L在零抖动下选择已记录细化增益最大的点；仅用于诊断", ...
    'experimentSource',fileread(mfilename('fullpath')+".m"));
directory = fileparts(resultFile);
save(fullfile(directory,"direction_only_control_"+cfg.profile+".mat"),'result');
writetable(result.records,fullfile(directory,"direction_only_control_"+cfg.profile+".csv"));
disp(result.records);
end
