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
            representedErrors = zeros(size(gridSizes));
            halfWidth = deg2rad(cfg.spectrum.desiredExclusionHalfWidthDeg);
            expectedIncludedSolidAngle = 2 * pi * cos(halfWidth);

            for index = 1:numel(gridSizes)
                gridSize = gridSizes(index);
                solidAngleGrid = uav.spectrum.upaIntegrationGrid(cfg, ...
                    gridSize, "solid_angle");
                azElGrid = uav.spectrum.upaIntegrationGrid(cfg, ...
                    gridSize, "az_el_solid_angle");
                testCase.verifyEqual(sum(solidAngleGrid.weights), 2 * pi, ...
                    'AbsTol', 2048 * eps(2 * pi));
                azElErrors(index) = abs(sum(azElGrid.weights) - 2 * pi);
                representedErrors(index) = abs( ...
                    solidAngleGrid.representedSolidAngle - ...
                    expectedIncludedSolidAngle);
            end

            testCase.verifyLessThan(azElErrors(2:end), azElErrors(1:end-1));
            testCase.verifyLessThanOrEqual(representedErrors, ...
                [0.08, 0.02, 0.01, 0.01]);
            testCase.verifyLessThan(representedErrors(end), representedErrors(1));
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
    end
end
