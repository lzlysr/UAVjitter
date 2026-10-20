classdef TestRecRootCauseControls < matlab.unittest.TestCase
    %TESTRECROOTCAUSECONTROLS 锁定原式迁移诊断的解析依据及反例。
    methods (TestMethodSetup)
        function addPaths(testCase)
            root = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root, 'src')));
        end
    end
    methods (Test)
        function refinedPointModelRepairsPopulationGeometryCase(testCase)
            % 总体协方差干预：X只是接口维度占位，不是有限64快拍性能结论。
            % 同一零失配输入下，验证新点源估计器消除了旧二维积分的大部分损失。
            cfg = fixtureConfig();
            candidate = cfg.experiments.recRootCause.candidateCovfit;
            cfg.spectrum.gridSize = candidate.gridSize;
            cfg.literature.guLeshem2012.upaDiagnosis.loadingRatio = candidate.loadingRatio;
            cfg.literature.guLeshem2012.upaDiagnosis.peakNmsRadiusBins = candidate.peakNmsRadiusBins;
            cfg.literature.guLeshem2012.upaDiagnosis.fitRegularization = candidate.fitRegularization;
            cfg.literature.guLeshem2012.upaDiagnosis.refineCandidates = true;
            cfg = uav.config.validate(cfg);
            A = cfg.derived.nominalSteering;
            truth = eye(9)+A(:,2:end)*diag(cfg.derived.interfererPowers)*A(:,2:end)';
            Rx = truth+cfg.derived.desiredPower*A(:,1)*A(:,1)';
            safe = struct('X',zeros(9,64),'Rx',Rx);
            [actual,d] = uav.incm.upaCovfit(safe,cfg);
            w = uav.beamforming.solve("mvdr",actual,A(:,1),cfg);
            value = 10*log10(cfg.derived.desiredPower*abs(w'*A(:,1))^2/real(w'*truth*w));
            testCase.verifyGreaterThan(value,19.3);
            testCase.verifyLessThan(norm(actual-truth,'fro')^2/norm(truth,'fro')^2,1e-4);
            testCase.verifyFalse(d.usedTruth);
        end
        function refinementKeepsRectangularArrayOrdering(testCase)
            % Nx不等于Ny时仍应使用y索引先变化的阵元顺序。
            cfg = fixtureConfig();
            cfg.array.Nx = 2;
            cfg.array.Ny = 3;
            cfg = uav.config.validate(cfg);
            A = cfg.derived.nominalSteering;
            R = eye(6) + A*diag(cfg.derived.sourcePowers)*A';
            nodes = cfg.derived.nominalDirectionCosines(:,2:end).' + [0.003,-0.004];
            [actual,d] = uav.spectrum.refineCaponNodes(nodes,R,false,cfg);
            testCase.verifySize(actual,[2,2]);
            testCase.verifyLessThan(max(abs(actual-cfg.derived.nominalDirectionCosines(:,2:end).'),[],'all'),0.01);
            testCase.verifyTrue(all(d.accepted));
        end
        function refinedCovfitDoesNotReadConfiguredInterfererDirections(testCase)
            % 将非期望源的配置派生字段置为无意义值，不改变观测与期望先验。
            % 若算法偷读真实干扰方向或数量，本测试将报错或改变输出。
            cfg = fixtureConfig();
            cfg.timeline.L = 64;
            cfg.literature.guLeshem2012.upaDiagnosis.refineCandidates = true;
            scene = uav.scenario.generate(cfg,uav.util.makeStream(cfg.random.baseSeed,11));
            observed = uav.scenario.generatePeriod(cfg,scene,1,uav.util.makeStream(cfg.random.baseSeed,12));
            safe = uav.experiment.algorithmObservation(observed);
            [first,d] = uav.incm.upaCovfit(safe,cfg);
            poisoned = cfg;
            poisoned.sources.interfererAzElDeg = [];
            poisoned.derived.J = 99;
            poisoned.derived.sourcePowers = [];
            poisoned.derived.nominalArrayUnitVectors(:,2:end) = NaN;
            poisoned.derived.nominalDirectionCosines(:,2:end) = NaN;
            poisoned.derived.nominalSteering(:,2:end) = NaN;
            second = uav.incm.upaCovfit(safe,poisoned);
            testCase.verifyEqual(first,second,'AbsTol',0);
            testCase.verifyTrue(d.candidateRefinement.enabled);
            testCase.verifyFalse(d.usedTruth);
        end
        function matchedWorldSourcesAlsoWorkWithNineSensorUla(testCase)
            % 同源功率、同世界方向及同阵元数；只改为y轴线阵及其可观测一维角域。
            % 线阵孔径与三维先验含义随之变化，不能将本对照称为纯维度效应。
            cfg = fixtureConfig();
            cfg.array.Nx = 1;
            cfg.array.Ny = 9;
            cfg = uav.config.validate(cfg);
            A = cfg.derived.nominalSteering;
            Ri = eye(9) + A(:,2:end)*diag(cfg.derived.interfererPowers)*A(:,2:end)';
            Rx = Ri + cfg.derived.desiredPower*A(:,1)*A(:,1)';
            step = pi/4001;
            theta = -pi/2 + step/2 + step*(0:4000);
            allowed = abs(theta-asin(cfg.derived.nominalDirectionCosines(2,1))) > deg2rad(25);
            G = uav.array.steeringUPA([zeros(1,sum(allowed)); sin(theta(allowed))], ...
                1,9,[0.5,0.5]);
            R = uav.incm.guLeshem2012Rec(Rx,G,step*ones(1,sum(allowed)));
            [w,d] = uav.beamforming.solve("mvdr",R,A(:,1),cfg);
            value = 10*log10(cfg.derived.desiredPower*abs(w'*A(:,1))^2/real(w'*Ri*w));
            fprintf('MATCHED_WORLD_ULA SINR=%.9f\n',value);
            testCase.verifyGreaterThan(value,15);
            testCase.verifyEqual(d.loading,0);
        end
        function observationOnlyRefinementRecoversOffGridInterferers(testCase)
            % 初始节点人为偏离真方向，细化函数只接收总协方差与先验角域。
            cfg = fixtureConfig();
            A = cfg.derived.nominalSteering;
            R = eye(9) + A * diag(cfg.derived.sourcePowers) * A';
            nodes = cfg.derived.nominalDirectionCosines(:, 2:end).' + [0.003, -0.004];
            [refined, d] = uav.spectrum.refineCaponNodes(nodes, R, false, cfg);
            testCase.verifyLessThan(max(abs(refined - ...
                cfg.derived.nominalDirectionCosines(:, 2:end).'), [], 'all'), 1e-3);
            testCase.verifyTrue(all(d.exitflags > 0));
            testCase.verifyFalse(d.usedTruth);
        end
        function singleSourceCaponMatchesAnalyticFormula(testCase)
            % 不经过样本生成，检查峰高、峰宽的解析公式。
            M = 9;
            power = 100;
            noise = 1;
            a = uav.array.steeringUPA([0; 0], 3, 3, [0.5, 0.5]);
            G = uav.array.steeringUPA([0, 0.2, -0.4; 0, -0.3, 0.1], ...
                3, 3, [0.5, 0.5]);
            R = noise * eye(M) + power * a * a';
            [~, diagnostics] = uav.incm.guLeshem2012Rec(R, G, ones(1, 3));
            similarity = abs(a' * G).^2 / M^2;
            expected = noise / M * (noise + M * power) ./ ...
                (noise + M * power * (1 - similarity));
            testCase.verifyEqual(diagnostics.caponPower, expected, 'RelTol', 1e-11);
            testCase.verifyEqual(diagnostics.caponPower(1), power + noise/M, 'AbsTol', 1e-10);
        end
        function pureNoiseIntegralIsNotAutomaticallyWhite(testCase)
            % 立体角下，平面阵对角方向的空间相关不必为零。
            cfg = fixtureConfig();
            grid = uav.spectrum.upaIntegrationGrid(cfg, 101, "solid_angle");
            [R, d] = uav.incm.guLeshem2012Rec(eye(9), grid.manifold, grid.weights);
            expected = (grid.manifold .* (grid.weights / 9)) * grid.manifold';
            testCase.verifyEqual(R, expected, 'AbsTol', 1e-11);
            testCase.verifyEqual(d.loading, 0);
            offDiagonal = R - diag(diag(R));
            testCase.verifyGreaterThan(norm(offDiagonal, 'fro') / norm(R, 'fro'), 0.1);
        end
        function exactPopulationStillShowsDefaultGeometryLoss(testCase)
            cfg = fixtureConfig();
            [value, upper, loading] = populationSinr(cfg);
            testCase.verifyLessThan(value, 11);
            testCase.verifyGreaterThan(upper, 19);
            testCase.verifyEqual(loading, 0);
        end
        function sameUpaSupportsSuccessfulOriginalRecCounterexample(testCase)
            cfg = fixtureConfig();
            cfg.sources.desiredAzElDeg = [0, -90];
            cfg.sources.interfererAzElDeg = [0, -acosd(2/3); 90, -acosd(2/3)];
            cfg = uav.config.validate(cfg);
            [value, upper, loading] = populationSinr(cfg);
            testCase.verifyGreaterThan(value, 18);
            testCase.verifyLessThan(upper - value, 1);
            testCase.verifyEqual(loading, 0);
        end
    end
end

function cfg = fixtureConfig()
%FIXTURECONFIG 人工固定几何与功率，防止用户日常配置改变诊断合同。
cfg = config('smoke');
cfg.array.Nx = 3;
cfg.array.Ny = 3;
cfg.array.spacingOverLambda = [0.5, 0.5];
cfg.attitude.nominalDeg = [0, 0, 0];
cfg.attitude.stdDeg = [0, 0, 0];
cfg.sources.desiredAzElDeg = [0, -35];
cfg.sources.interfererAzElDeg = [40, -25; -55, -20];
cfg.sources.interfererInrDb = [30, 30];
cfg.signal.noisePower = 1;
cfg.signal.snrDb = 10;
cfg.spectrum.desiredExclusionHalfWidthDeg = 25;
cfg = uav.config.validate(cfg);
end

function [value, upper, loading] = populationSinr(cfg)
%POPULATIONSINR 仅作理论诊断，不调用待测试的实验包装或样本生成器。
A = cfg.derived.nominalSteering;
Ri = eye(9) + A(:, 2:end) * diag(cfg.derived.interfererPowers) * A(:, 2:end)';
Rx = Ri + cfg.derived.desiredPower * A(:, 1) * A(:, 1)';
grid = uav.spectrum.upaIntegrationGrid(cfg, 201, "solid_angle");
mask = grid.includedMask;
R = uav.incm.guLeshem2012Rec(Rx, grid.manifold(:, mask), grid.weights(mask));
[w, diagnostics] = uav.beamforming.solve("mvdr", R, A(:, 1), cfg);
wo = uav.beamforming.solve("mvdr", Ri, A(:, 1), cfg);
value = 10 * log10(cfg.derived.desiredPower * abs(w' * A(:, 1))^2 / real(w' * Ri * w));
upper = 10 * log10(cfg.derived.desiredPower * abs(wo' * A(:, 1))^2 / real(wo' * Ri * wo));
loading = diagnostics.loading;
end
