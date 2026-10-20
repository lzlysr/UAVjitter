function results = runCovfitPopulationDiagnosis(profile)
%RUNCOVFITPOPULATIONDIAGNOSIS 检查既有CovFit能否在总体极限恢复有效INCM。
%   只研究已实现估计器，参数全部交叉而非只报告有利点。
%   真实协方差仅用于统一输入干预；候选生成及功率估计不读取真实方向。
arguments
    profile (1,1) string = "smoke"
end
cfg = config(profile);
study = cfg.experiments.recRootCause;
cfg.attitude.stdDeg = zeros(size(cfg.attitude.stdDeg));
% 大L只用于选择既有噪声估计分支；Rx在下方替换为精确总体输入。
cfg.timeline.L = max(study.snapshotCounts);
cfg = uav.config.validate(cfg);
stream = uav.util.makeStream(cfg.random.baseSeed, 1);
scenario = uav.scenario.generate(cfg, stream);
obs = uav.scenario.generatePeriod(cfg, scenario, 1, stream);
A = obs.truth.trueASVs;
Ri = obs.truth.trueIncm;
obs.Rx = Ri + cfg.derived.desiredPower * A(:, 1) * A(:, 1)';
safe = uav.experiment.algorithmObservation(obs);
records = struct([]);
timer = tic;
for gridSize = study.covfitGridSizes
    local = cfg;
    local.spectrum.gridSize = gridSize;
    local = uav.config.validate(local);
    for loading = study.covfitLoadingRatios
        for radius = study.covfitNmsRadii
            for regularization = study.covfitRegularizations
                local.literature.guLeshem2012.upaDiagnosis.loadingRatio = loading;
                local.literature.guLeshem2012.upaDiagnosis.peakNmsRadiusBins = radius;
                local.literature.guLeshem2012.upaDiagnosis.fitRegularization = regularization;
                [Rhat, d] = uav.incm.upaCovfit(safe, local);
                [w, wd] = uav.beamforming.solve("mvdr", Rhat, ...
                    cfg.derived.nominalSteeringReference, local);
                record = struct('gridSize', gridSize, 'loadingRatio', loading, ...
                    'nmsRadius', radius, 'regularization', regularization, ...
                    'sinrDb', 10*log10(cfg.derived.desiredPower * abs(w'*A(:,1))^2 / real(w'*Ri*w)), ...
                    'incmNmse', norm(Rhat - Ri, 'fro')^2 / norm(Ri, 'fro')^2, ...
                    'fitResidual', d.fitResidual, 'candidateCount', d.interfererCandidateCount, ...
                    'interferencePowerSum', sum(d.interfererFittedPowers), ...
                    'noiseEstimate', d.noisePower, 'weightLoading', wd.loading, ...
                    'usedTruth', d.usedTruth);
                records = [records; record]; %#ok<AGROW>
            end
        end
    end
end
root = fileparts(fileparts(mfilename('fullpath')));
outputDirectory = fullfile(root, study.outputDirectory);
if ~isfolder(outputDirectory), mkdir(outputDirectory); end
results = struct('cfg', cfg, 'records', struct2table(records), ...
    'metadata', uav.experiment.metadata(cfg, "covfit-population-diagnosis", ...
        struct('parameterPoints', numel(records)), toc(timer)), ...
    'experimentSource', fileread(mfilename('fullpath') + ".m"), ...
    'configurationSource', fileread(fullfile(root, 'config.m')), ...
    'inputCovariance', obs.Rx, 'diagnosticOnlyTruth', true);
save(fullfile(outputDirectory, 'covfit_population.mat'), 'results', '-v7.3');
writetable(results.records, fullfile(outputDirectory, 'covfit_population.csv'));
disp(sortrows(results.records, 'sinrDb', 'descend'));
end
