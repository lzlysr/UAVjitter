function result = runRawRecNarrowControl(resultFile)
%RUNRAWRECNARROWCONTROL 对原式只改变角域，重放完整验证的同一批快拍。
%   不改变积分测度、网格、加载或权值参考；L<M仍标记为不适用。
arguments
    resultFile (1,1) string
end
saved = load(resultFile,'results');
r = saved.results;
cfg = r.cfg;
s = cfg.experiments.recRootCause;
cfg.spectrum.desiredExclusionHalfWidthDeg = s.narrowSectorHalfWidthDeg;
grid = uav.spectrum.upaIntegrationGrid(cfg,s.diagnosticGridSize,"solid_angle");
mask = grid.includedMask;
G = grid.manifold(:,mask);
weights = grid.weights(mask);
sinrLinear = nan(s.validationTrials,s.validationPeriods, ...
    numel(s.validationJitterDeg),numel(s.validationSnapshotCounts));
weightLoading = sinrLinear;
timer = tic;
for trial = 1:s.validationTrials
    for point = 1:numel(s.validationJitterDeg)
        local = cfg;
        local.attitude.stdDeg(:) = s.validationJitterDeg(point);
        scene = uav.scenario.generate(local,uav.util.makeStream(cfg.random.baseSeed,s.validationSeedOffset+trial));
        for k = 1:s.validationPeriods
            obs = uav.scenario.generatePeriod(local,scene,k,uav.util.makeStream( ...
                cfg.random.baseSeed,s.validationSeedOffset+s.validationTrials+(trial-1)*s.validationPeriods+k));
            for snapshotIndex = 1:numel(s.validationSnapshotCounts)
                L = s.validationSnapshotCounts(snapshotIndex);
                if L<cfg.derived.M, continue; end
                Rx = obs.X(:,1:L)*obs.X(:,1:L)'/L;
                Rh = uav.incm.guLeshem2012Rec(Rx,G,weights);
                [w,d] = uav.beamforming.solve("mvdr",Rh,cfg.derived.nominalSteeringReference,cfg);
                m = uav.metrics.evaluate(w,obs.truth,cfg);
                sinrLinear(trial,k,point,snapshotIndex) = m.sinrLinear;
                weightLoading(trial,k,point,snapshotIndex) = d.loading;
            end
        end
    end
    fprintf('RAW_NARROW_CONTROL trial=%d/%d elapsed=%.1fs\n',trial,s.validationTrials,toc(timer));
end
rows = struct([]);
originalIndex = r.methodNames=="gu_leshem2012_rec";
targetIndex = r.methodNames=="refined_narrow_sector";
for snapshotIndex = 1:numel(s.validationSnapshotCounts)
    for point = 1:numel(s.validationJitterDeg)
        raw = sinrLinear(:,:,point,snapshotIndex);
        original = r.sinrLinear(:,:,point,snapshotIndex,originalIndex);
        candidate = r.sinrLinear(:,:,point,snapshotIndex,targetIndex);
        row = struct('L',s.validationSnapshotCounts(snapshotIndex),'jitterDeg',s.validationJitterDeg(point), ...
            'rawRec25MedianDb',10*log10(median(mean(original,2))), ...
            'rawRec7MedianDb',10*log10(median(mean(raw,2))), ...
            'refinedCovfit7MedianDb',10*log10(median(mean(candidate,2))), ...
            'rawRec7P05Db',10*log10(prctile(raw(:),5)));
        rows = [rows;row]; %#ok<AGROW>
    end
end
result = struct('cfg',cfg,'sourceValidation',resultFile,'sinrLinear',sinrLinear, ...
    'weightLoading',weightLoading,'summary',struct2table(rows),'runtimeSec',toc(timer), ...
    'experimentSource',fileread(mfilename('fullpath')+".m"));
directory = fileparts(resultFile);
save(fullfile(directory,"raw_rec_narrow_control_"+cfg.profile+".mat"),'result','-v7.3');
writetable(result.summary,fullfile(directory,"raw_rec_narrow_control_"+cfg.profile+".csv"));
disp(result.summary);
end
