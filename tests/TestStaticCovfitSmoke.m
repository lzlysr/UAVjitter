classdef TestStaticCovfitSmoke < matlab.unittest.TestCase
    %TESTSTATICCOVFITSMOKE 验证零抖动smoke场景中的经典静态INCM替换。

    methods (Test)
        function classicStaticMethodUsesUpaCovfitAtZeroJitter(testCase)
            % 若经典方法仍调用旧的均匀Capon积分，本测试会因估计器标识缺失而失败。
            cfg = config("smoke");
            cfg.timeline.L = 8;
            cfg.attitude.stdDeg = [0, 0, 0];
            cfg = uav.config.validate(cfg);

            scenario = uav.scenario.generate(cfg, ...
                uav.util.makeStream(cfg.random.baseSeed, 901));
            obs = uav.scenario.generatePeriod(cfg, scenario, 1, ...
                uav.util.makeStream(cfg.random.baseSeed, 902));
            [record, ~] = uav.experiment.runClassicMethods( ...
                obs, [], cfg, "static_incm_mvdr");

            reconstruction = record.diagnostics{1}.reconstruction;
            testCase.verifyEqual(reconstruction.estimator, "upa-covfit");
            testCase.verifyFalse(reconstruction.usedTruth);
            testCase.verifyTrue(all(isfinite(record.weights), 'all'));
            testCase.verifyTrue(isfinite(record.sinrLinear));
            testCase.verifyGreaterThan(record.sinrLinear, 0);
            testCase.verifyGreaterThanOrEqual(record.incmMinimumEigenvalue, -1e-10);
        end
    end
end
