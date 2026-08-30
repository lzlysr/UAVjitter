classdef TestTask8UpaIncmDiagnosis < matlab.unittest.TestCase
    %TESTTASK8UPAINCMDIAGNOSIS 验证UPA诊断配置和三类积分网格合同。

    methods (TestMethodSetup)
        function addRuntimePaths(testCase)
            % 测试直接调用工程包函数和公开配置入口。
            rootDir = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(rootDir));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(rootDir, 'src')));
        end
    end

    methods (Test)
        % 防止协方差拟合把两个候选源的独立功率混合或读取任何真值。
        function nonnegativeFitRecoversSyntheticSourcePowers(testCase)
            cfg = config("smoke");
            xi = [0.2, -0.4; -0.3, 0.1];
            A = uav.array.steeringUPA(xi, cfg.array.Nx, cfg.array.Ny, ...
                cfg.array.spacingOverLambda);
            expected = [20; 7];
            noise = 1;
            Rx = A * diag(expected) * A' + noise * eye(cfg.derived.M);

            [actual, diagnostics] = uav.incm.fitNonnegativePowers( ...
                Rx, A, noise, 0);

            testCase.verifyEqual(actual, expected, 'RelTol', 1e-9);
            testCase.verifyLessThan(diagnostics.relativeResidual, 1e-10);
            testCase.verifyEqual(diagnostics.dictionaryRank, 2);
            testCase.verifyFalse(diagnostics.usedTruth);
        end

        % 防止严格下三角未按sqrt(2)缩放，从而改变Frobenius正则化权衡。
        function nonnegativeFitPreservesFrobeniusWeightForOneCandidate(testCase)
            cfg = config("smoke");
            A = uav.array.steeringUPA([0.2; -0.3], cfg.array.Nx, ...
                cfg.array.Ny, cfg.array.spacingOverLambda);
            expectedSourcePower = 20;
            lambdaFit = 9;
            Rx = expectedSourcePower * (A * A') + eye(cfg.derived.M);
            expected = 18;

            actual = uav.incm.fitNonnegativePowers(Rx, A, 1, lambdaFit);

            testCase.verifyEqual(actual, expected, 'AbsTol', 1e-10);
        end

        % 防止非法协方差、空候选或非法标量悄然进入NNLS求解器。
        function nonnegativeFitRejectsInvalidInputs(testCase)
            validCovariance = eye(2);
            validManifold = [1; 1];

            testCase.verifyError(@() uav.incm.fitNonnegativePowers( ...
                ones(2, 3), validManifold, 0, 0), ...
                'uav:incm:InvalidCovariance');
            testCase.verifyError(@() uav.incm.fitNonnegativePowers( ...
                validCovariance, zeros(2, 0), 0, 0), ...
                'uav:incm:InvalidCandidateManifold');
            testCase.verifyError(@() uav.incm.fitNonnegativePowers( ...
                validCovariance, validManifold, -1, 0), ...
                'uav:incm:InvalidNoisePower');
            testCase.verifyError(@() uav.incm.fitNonnegativePowers( ...
                validCovariance, validManifold, 0, -1), ...
                'uav:incm:InvalidFitRegularization');
            testCase.verifyError(@() uav.incm.fitNonnegativePowers( ...
                [1, NaN; 0, 1], validManifold, 0, 0), ...
                'uav:incm:InvalidCovariance');
        end

        % 防止诊断配置遗漏低快拍条件或任一指定求积测度。
        function configurationContainsBothSnapshotCounts(testCase)
            cfg = config("smoke");
            d = cfg.literature.guLeshem2012.upaDiagnosis;
            testCase.verifyEqual(d.snapshotCounts, [8, 64]);
            testCase.verifyEqual(d.gridSizes, [21, 51, 101, 201]);
            testCase.verifyEqual(d.measureNames, ...
                ["direction_cosine", "solid_angle", "az_el_solid_angle"]);
        end

        % 防止积分网格返回非有限或非正权重，或把期望排除区重新纳入积分。
        function integrationGridsUsePositiveFiniteWeights(testCase)
            cfg = config("smoke");
            for name = cfg.literature.guLeshem2012.upaDiagnosis.measureNames
                grid = uav.spectrum.upaIntegrationGrid(cfg, 51, name);
                testCase.verifySize(grid.manifold, ...
                    [cfg.derived.M, size(grid.nodes, 1)]);
                testCase.verifyTrue(all(isfinite(grid.weights)));
                testCase.verifyGreaterThan(min(grid.weights), 0);
                testCase.verifyFalse(any(grid.includedMask & grid.desiredMask));
            end
        end

        % 防止solid_angle测度在方向余弦圆盘边界出现雅可比奇异或错误半球面积。
        function solidAngleMeasuresRepresentTheLowerHemisphere(testCase)
            cfg = config("smoke");
            gridSizes = [21, 51, 101, 201];
            azElErrors = zeros(size(gridSizes));
            representedErrors = zeros(3, numel(gridSizes));
            halfWidth = deg2rad(cfg.spectrum.desiredExclusionHalfWidthDeg);
            expectedIncludedSolidAngle = 2 * pi * cos(halfWidth);

            for index = 1:numel(gridSizes)
                gridSize = gridSizes(index);
                directionCosineGrid = uav.spectrum.upaIntegrationGrid(cfg, ...
                    gridSize, "direction_cosine");
                solidAngleGrid = uav.spectrum.upaIntegrationGrid(cfg, ...
                    gridSize, "solid_angle");
                azElGrid = uav.spectrum.upaIntegrationGrid(cfg, ...
                    gridSize, "az_el_solid_angle");
                testCase.verifyEqual(sum(solidAngleGrid.weights), 2 * pi, ...
                    'AbsTol', 2048 * eps(2 * pi));
                azElErrors(index) = abs(sum(azElGrid.weights) - 2 * pi);
                grids = [directionCosineGrid, solidAngleGrid, azElGrid];
                for measureIndex = 1:numel(grids)
                    representedErrors(measureIndex, index) = abs( ...
                        grids(measureIndex).representedSolidAngle - ...
                        expectedIncludedSolidAngle);
                end
            end

            testCase.verifyLessThan(azElErrors(2:end), azElErrors(1:end-1));
            testCase.verifyLessThanOrEqual(representedErrors(1:2, :), ...
                repmat([0.08, 0.02, 0.01, 0.01], 2, 1));
            testCase.verifyLessThanOrEqual(representedErrors(3, :), ...
                [0.08, 0.02, 0.01, 0.01]);
            testCase.verifyLessThan(representedErrors(1, end), ...
                representedErrors(1, 1));
        end

        % 防止诊断配置把向量误当单值参数，或放宽低快拍和正加载合同。
        function diagnosisConfigurationRejectsInvalidScalarContracts(testCase)
            invalidCases = { ...
                'snapshotCounts', [9, 64]; ...
                'loadingRatio', 0; ...
                'loadingRatio', [1e-2, 2e-2]; ...
                'fitRegularization', -1e-4; ...
                'fitRegularization', [0, 1e-4]; ...
                'peakThresholdRatio', 1.1; ...
                'peakThresholdRatio', [0.1, 0.2]; ...
                'peakNmsRadiusBins', [1, 2]; ...
                'maximumCandidateCount', 0; ...
                'maximumCandidateCount', [8, 9]; ...
                'monteCarloTrials', 0; ...
                'monteCarloTrials', [1, 2]; ...
                'periods', 0; ...
                'periods', [1, 2]; ...
                'experimentIndex', 0; ...
                'experimentIndex', [9, 10]; ...
                'selectedGridSize', [101, 201]; ...
                'selectedMeasureName', ["solid_angle", "direction_cosine"]};
            for index = 1:size(invalidCases, 1)
                cfg = config("smoke");
                field = invalidCases{index, 1};
                cfg.literature.guLeshem2012.upaDiagnosis.(field) = ...
                    invalidCases{index, 2};
                testCase.verifyError(@() uav.config.validate(cfg), ...
                    'uav:config:InvalidExperimentParameter');
            end

            cfg = config("smoke");
            cfg.literature.guLeshem2012.upaDiagnosis.fitRegularization = 0;
            testCase.verifyWarningFree(@() uav.config.validate(cfg));
        end

        % 防止紧凑诊断被迫计算未选择的测度，同时拒绝空集和重复测度。
        function diagnosisConfigurationAllowsSupportedMeasureSubset(testCase)
            cfg = config("smoke");
            cfg.literature.guLeshem2012.upaDiagnosis.measureNames = ...
                ["direction_cosine", "solid_angle"];
            cfg.literature.guLeshem2012.upaDiagnosis.selectedMeasureName = ...
                "solid_angle";
            validated = uav.config.validate(cfg);
            testCase.verifyEqual( ...
                validated.literature.guLeshem2012.upaDiagnosis.measureNames, ...
                ["direction_cosine", "solid_angle"]);

            invalid = cfg;
            invalid.literature.guLeshem2012.upaDiagnosis.measureNames = ...
                ["solid_angle", "solid_angle"];
            testCase.verifyError(@() uav.config.validate(invalid), ...
                'uav:config:InvalidExperimentParameter');
        end

        % 防止欠定快拍仍误用过时配置L，或在重构中泄漏真实方向与真实INCM。
        function upaCovfitUsesObservedUnderdeterminedSnapshotsWithoutTruth(testCase)
            cfg = config("smoke");
            cfg.timeline.L = 64;
            cfg = uav.config.validate(cfg);
            scenario = uav.scenario.generate(cfg, uav.util.makeStream(211, 1));
            obs = uav.scenario.generatePeriod(cfg, scenario, 1, ...
                uav.util.makeStream(211, 2));
            obs.X = obs.X(:, 1:8);
            obs.pilot = obs.pilot(1:8);
            obs.Rx = obs.X * obs.X' / size(obs.X, 2);
            safe = uav.experiment.algorithmObservation(obs);
            originalLoading = cfg.spectrum.diagonalLoading;

            [RinHat, diagnostics] = uav.incm.upaCovfit(safe, cfg);

            testCase.verifyFalse(isfield(safe, 'truth'));
            testCase.verifyTrue(all(isfinite(RinHat), 'all'));
            testCase.verifyLessThan(norm(RinHat - RinHat', 'fro'), 1e-10);
            testCase.verifyGreaterThan(min(real(eig(RinHat))), 0);
            testCase.verifyEqual(diagnostics.noiseSource, ...
                "configured-underdetermined");
            testCase.verifyEqual(diagnostics.effectiveSnapshotCount, 8);
            testCase.verifyEqual(diagnostics.noisePower, cfg.signal.noisePower);
            testCase.verifyEqual(cfg.spectrum.diagonalLoading, originalLoading);
            testCase.verifyFalse(diagnostics.usedTruth);
        end

        % 防止名义期望与期望区谱峰重复进入字典，或干扰候选超过阵列秩上限。
        function upaCovfitDeduplicatesCandidatesAndRecordsAuditableFit(testCase)
            cfg = config("smoke");
            scenario = uav.scenario.generate(cfg, uav.util.makeStream(223, 1));
            obs = uav.scenario.generatePeriod(cfg, scenario, 1, ...
                uav.util.makeStream(223, 2));
            safe = uav.experiment.algorithmObservation(obs);

            [RinHat, diagnostics] = uav.incm.upaCovfit(safe, cfg);

            candidateNodes = [diagnostics.desiredCandidateNodes; ...
                diagnostics.interfererCandidateNodes];
            tolerance = cfg.spectrum.quadratureDedupTolerance;
            testCase.verifyEqual(size(unique(round(candidateNodes / tolerance), ...
                'rows'), 1), size(candidateNodes, 1));
            testCase.verifyLessThanOrEqual(size(diagnostics.interfererCandidateNodes, 1), ...
                min(cfg.literature.guLeshem2012.upaDiagnosis.maximumCandidateCount, ...
                cfg.derived.M - 1));
            testCase.verifySize(diagnostics.fittedPowers, ...
                [size(candidateNodes, 1), 1]);
            testCase.verifyTrue(isfinite(diagnostics.fitResidual));
            testCase.verifyTrue(isfinite(diagnostics.conditionNumber));
            testCase.verifyTrue(all(isfinite(RinHat), 'all'));
        end

        % 防止M=1时无峰回退绕过min(maximumCandidateCount,M-1)的硬上限。
        function upaCovfitKeepsNoInterfererCandidatesWhenArrayHasOneSensor(testCase)
            cfg = config("smoke");
            cfg.array.Nx = 1;
            cfg.array.Ny = 1;
            cfg = uav.config.validate(cfg);
            scenario = uav.scenario.generate(cfg, uav.util.makeStream(227, 1));
            obs = uav.scenario.generatePeriod(cfg, scenario, 1, ...
                uav.util.makeStream(227, 2));
            safe = uav.experiment.algorithmObservation(obs);

            [RinHat, diagnostics] = uav.incm.upaCovfit(safe, cfg);

            testCase.verifyGreaterThanOrEqual(cfg.derived.J, cfg.derived.M);
            testCase.verifyEqual(diagnostics.interfererCandidateCap, 0);
            testCase.verifyEqual(diagnostics.interfererCandidateCount, 0);
            testCase.verifyEqual(RinHat, diagnostics.noisePower * eye(1), ...
                'AbsTol', 1e-12);
            testCase.verifyTrue(all(isfinite(RinHat), 'all'));
            testCase.verifyGreaterThan(min(real(eig(RinHat))), 0);
        end

        % 防止名义期望与精确谱网格峰重合时仍形成重复的期望拟合字典列。
        function upaCovfitDeduplicatesNominalNodeAndDesiredGridPeak(testCase)
            cfg = config("smoke");
            cfg.sources.desiredAzElDeg = [0, -90];
            cfg = uav.config.validate(cfg);
            nominal = cfg.derived.nominalSteeringReference;
            snapshotCount = cfg.timeline.L;
            X = sqrt(30) * nominal * ones(1, snapshotCount) + ...
                sqrt(cfg.signal.noisePower) * eye(cfg.derived.M, snapshotCount);
            safe = struct('X', X, 'Rx', X * X' / snapshotCount);

            [~, diagnostics] = uav.incm.upaCovfit(safe, cfg);

            testCase.verifyEqual(size(diagnostics.desiredCandidateNodes, 1), 1);
            testCase.verifyEqual(diagnostics.desiredCandidateFallback, ...
                "deduplicated-with-nominal");
        end

        % 防止短快拍不是长快拍的逐元素前缀，或遗漏pilot和真值干扰符号。
        function shortObservationIsExactPrefixAndRecomputesCovariance(testCase)
            [~, obs64] = testCase.zeroJitterObservation(64);

            obs8 = uav.scenario.truncatePeriodObservation(obs64, 8);

            testCase.verifyEqual(obs8.X, obs64.X(:, 1:8));
            testCase.verifyEqual(obs8.pilot, obs64.pilot(1:8));
            testCase.verifyEqual(obs8.truth.interfererSymbols, ...
                obs64.truth.interfererSymbols(:, 1:8));
            testCase.verifyEqual(obs8.Rx, obs8.X * obs8.X' / 8, ...
                'AbsTol', 1e-13);
            testCase.verifyEqual(obs64.Rx, obs64.X * obs64.X' / 64, ...
                'AbsTol', 1e-13);
            testCase.verifyEqual(obs8.truth.trueIncm, obs64.truth.trueIncm);
        end

        % 防止L<M时短快拍原始Rec仍调用核心函数或伪装为可用结果。
        function shortObservationMarksOriginalRecUnavailableWithoutCoreCall(testCase)
            [cfg64, obs64] = testCase.zeroJitterObservation(64);
            obs8 = uav.scenario.truncatePeriodObservation(obs64, 8);
            cfg8 = cfg64;
            cfg8.timeline.L = 8;
            cfg8 = uav.config.validate(cfg8);
            grid = uav.spectrum.upaIntegrationGrid(cfg8, 21, "solid_angle");

            record = uav.experiment.runUpaDiagnosisMethods(obs8, cfg8, ...
                ["oracle", "gu_leshem2012_rec", "upa_covfit_incm"], grid);

            testCase.verifyTrue(isnan(record.sinrLinear(2)));
            testCase.verifyTrue(all(isnan(record.weights(:, 2))));
            testCase.verifyTrue(isnan(record.incmNmse(2)));
            testCase.verifyEqual(record.status(2), "singular-not-applicable");
            testCase.verifyEqual(record.diagnostics{2}.coreInvoked, false);
            testCase.verifyEqual(record.reasonCode(2), "L-less-than-M");
            testCase.verifyThat(record.reason(2), ...
                matlab.unittest.constraints.ContainsSubstring("快拍数"));
            testCase.verifyTrue(isfinite(record.sinrLinear(3)));
            testCase.verifyFalse(record.diagnostics{3}.usedTruth);
        end

        % 防止执行器擅自重排方法，或把未选择的方法插入结果。
        function upaDiagnosisRunnerPreservesArbitraryMethodSubsetAndOrder(testCase)
            [cfg, obs] = testCase.zeroJitterObservation(64);
            grid = uav.spectrum.upaIntegrationGrid(cfg, 21, "solid_angle");
            methods = ["wco", "oracle", "esp"];

            record = uav.experiment.runUpaDiagnosisMethods(obs, cfg, methods, grid);

            testCase.verifyEqual(record.methodNames, methods);
            testCase.verifySize(record.weights, [cfg.derived.M, numel(methods)]);
            testCase.verifyTrue(all(isfinite(record.sinrLinear)));
            testCase.verifyEqual(record.status, repmat("ok", 1, numel(methods)));
        end

        % 防止非Oracle方法在执行器中获得truth，或在指标之前读取真实量。
        function upaDiagnosisRunnerAuditsTruthIsolation(testCase)
            [cfg, obs] = testCase.zeroJitterObservation(64);
            grid = uav.spectrum.upaIntegrationGrid(cfg, 21, "solid_angle");
            methods = ["esp", "wco", "gu_leshem2012_rec", "upa_covfit_incm"];

            record = uav.experiment.runUpaDiagnosisMethods(obs, cfg, methods, grid);

            for methodIndex = 1:numel(methods)
                diagnostics = record.diagnostics{methodIndex};
                testCase.verifyFalse(diagnostics.algorithmInputHasTruth);
                testCase.verifyFalse(diagnostics.usedTruth);
            end
            testCase.verifyTrue(all(record.snapshotSignature == ...
                record.snapshotSignature(1, :), 'all'));
        end

        % 防止实验9遗漏配对维度、独立WCO诊断或任一种事务导出文件。
        function compactExperimentNineExportsPairedDiagnosisBundle(testCase)
            cfg = compactDiagnosisConfig(testCase);

            result = exp09_upa_incm_diagnosis(cfg);

            testCase.verifyEqual(result.methodNames, ...
                ["oracle", "esp", "wco", "gu_leshem2012_rec", ...
                "upa_covfit_incm"]);
            testCase.verifyEqual(result.snapshotCounts, [8, 64]);
            testCase.verifySize(result.raw.zeroJitter.sinrLinear, ...
                [1, 2, 2, 2, 5]);
            testCase.verifySize(result.raw.zeroJitter.status, ...
                [1, 2, 2, 2, 5]);
            testCase.verifySize(result.raw.jitterSweep.sinrLinear, ...
                [1, 1, 3, 2, 5]);
            testCase.verifySize( ...
                result.raw.jitterSweep.upaCovfitWcoSinrLinear, [1, 1, 3, 2]);
            testCase.verifyTrue(all( ...
                result.raw.zeroJitter.snapshotPrefixMatched, 'all'));
            testCase.verifyTrue(all( ...
                result.raw.jitterSweep.snapshotPrefixMatched, 'all'));
            testCase.verifyTrue(all(isnan( ...
                result.raw.zeroJitter.originalRec.traceRatio(:, 1, :, :)), ...
                'all'));
            testCase.verifyTrue(all(result.raw.zeroJitter.status( ...
                :, 1, :, :, 4) == "singular-not-applicable", 'all'));
            testCase.verifyTrue(all(isnan(result.raw.zeroJitter.sinrLinear( ...
                :, 1, :, :, 4)), 'all'));
            testCase.verifyTrue(all(isfinite( ...
                result.raw.zeroJitter.upaCovfit.traceRatio), 'all'));
            testCase.verifyEqual(result.raw.zeroJitter.sinrLinear(:, :, :, 1, 5), ...
                result.raw.zeroJitter.sinrLinear(:, :, :, 2, 5), 'AbsTol', 0);
            testCase.verifyEqual( ...
                result.summary.zeroJitter.sinr.aggregationDomain, "linear");
            testCase.verifyFalse(any(uav.experiment.entryPoints() == ...
                "exp09_upa_incm_diagnosis"));
            artifactPaths = string(struct2cell( ...
                rmfield(result.artifacts, 'status')));
            testCase.verifyTrue(all(isfile(artifactPaths)));
            testCase.verifyEqual(result.artifacts.status, "complete");
            testCase.verifyTrue(isfield(result, 'cfg'));
            testCase.verifyTrue(isfield(result, 'metadata'));
        end
    end

    methods (Access = private)
        function [cfg, obs] = zeroJitterObservation(testCase, snapshotCount)
            %ZEROJITTEROBSERVATION 构造同一零姿态抖动配对快拍的测试输入。
            cfg = config("smoke");
            cfg.timeline.L = snapshotCount;
            cfg.attitude.stdDeg = [0, 0, 0];
            cfg = uav.config.validate(cfg);
            scenario = uav.scenario.generate(cfg, uav.util.makeStream(251, 1));
            testCase.verifyEqual(scenario.attitudeOffsetRad, ...
                zeros(size(scenario.attitudeOffsetRad)), 'AbsTol', 1e-15);
            obs = uav.scenario.generatePeriod(cfg, scenario, 1, ...
                uav.util.makeStream(251, 2));
        end
    end
end


function cfg = compactDiagnosisConfig(testCase)
%COMPACTDIAGNOSISCONFIG 用隔离目录运行完整但最小的实验9链路。
cfg = config("smoke");
diagnosis = cfg.literature.guLeshem2012.upaDiagnosis;
diagnosis.monteCarloTrials = 1;
diagnosis.periods = 1;
diagnosis.gridSizes = [21, 31];
diagnosis.measureNames = ["direction_cosine", "solid_angle"];
diagnosis.selectedGridSize = 21;
diagnosis.selectedMeasureName = "solid_angle";
cfg.literature.guLeshem2012.upaDiagnosis = diagnosis;
cfg.experiments.bootstrapReplicates = 4;
cfg.plot.visible = false;
cfg.plot.closeAfterExport = true;
cfg.output.overwrite = true;
outputRoot = tempname;
mkdir(outputRoot);
testCase.addTeardown(@() cleanupDiagnosisDirectory(outputRoot));
cfg.output.rootDir = outputRoot;
cfg.output.rawDir = fullfile(outputRoot, 'raw');
cfg.output.tableDir = fullfile(outputRoot, 'tables');
cfg.output.figureDir = fullfile(outputRoot, 'figures');
cfg.output.logDir = fullfile(outputRoot, 'logs');
cfg.output.cacheDir = fullfile(outputRoot, 'cache');
cfg = uav.config.validate(cfg);
end

function cleanupDiagnosisDirectory(pathValue)
%CLEANUPDIAGNOSISDIRECTORY 清理本测试独占的可再生成导出文件。
if isfolder(pathValue)
    rmdir(pathValue, 's');
end
end
