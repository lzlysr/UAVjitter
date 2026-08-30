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
    end
end
