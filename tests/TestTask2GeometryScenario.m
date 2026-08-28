classdef TestTask2GeometryScenario < matlab.unittest.TestCase
    %TESTTASK2GEOMETRYSCENARIO 第二阶段：姿态、坐标、阵列和快拍模型测试。
    % 重点锁定旋转方向、传感器顺序、稳态姿态统计以及理论协方差，防止
    % “代码能运行但坐标系相反”这类难以从最终SINR发现的模型错误。
    methods (TestMethodSetup)
        % 为每个测试建立隔离的工程搜索路径。
        function addSourcePath(testCase)
            rootDir = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(rootDir));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(rootDir, 'src')));
        end
    end

    methods (Test)
        % R必须满足R'*R=I且det(R)=1，排除反射或非正交变换。
        function rotationZYXIsProperOrthogonal(testCase)
            R = uav.geometry.rotationZYX([pi / 6; -pi / 5; pi / 4]);

            testCase.verifyEqual(R' * R, eye(3), 'AbsTol', 1e-12);
            testCase.verifyEqual(det(R), 1, 'AbsTol', 1e-12);
        end

        % 零姿态时世界方向、阵列方向和方位/高度角应保持恒等映射。
        function zeroAttitudeKeepsWorldDirectionsAndCoordinates(testCase)
            uWorld = [1, 0; 0, 1; 0, 0];
            uArray = uav.geometry.worldToArray(eye(3), uWorld);
            [xi, azElRad] = uav.geometry.directionCoordinates(uArray);

            testCase.verifyEqual(uArray, uWorld, 'AbsTol', 1e-12);
            testCase.verifyEqual(xi, [1, 0; 0, 1], 'AbsTol', 1e-12);
            testCase.verifyEqual(azElRad, [0, pi / 2; 0, 0], 'AbsTol', 1e-12);
        end

        % 用手工可计算的非零ZYX旋转检查乘法顺序和坐标变换方向。
        function nonzeroRotationAndCoordinateConventionMatchesLiteralZYXCase(testCase)
            etaRad = [pi / 6; -pi / 4; pi / 3];
            R = uav.geometry.rotationZYX(etaRad);
            uWorld = [0; 0; -1];
            uArray = uav.geometry.worldToArray(R, uWorld);
            [xi, azElRad] = uav.geometry.directionCoordinates(uArray);
            expectedR = [ ...
                0.353553390593274, -0.926776695296637, 0.126826484044322; ...
                0.612372435695795, 0.126826484044322, -0.780330085889911; ...
                0.707106781186548, 0.353553390593274, 0.612372435695795];

            testCase.verifyEqual(R, expectedR, 'AbsTol', 1e-12);
            testCase.verifyEqual(uArray, [-0.707106781186548; ...
                -0.353553390593274; -0.612372435695795], 'AbsTol', 1e-12);
            testCase.verifyEqual(xi, [-0.707106781186548; -0.353553390593274], ...
                'AbsTol', 1e-12);
            testCase.verifyEqual(azElRad, [-2.677945044588987; -0.659058035826409], ...
                'AbsTol', 1e-12);
        end

        % 非单位方向、非法旋转和可见圆盘外方向必须被物理输入校验拒绝。
        function geometryAndSteeringRejectInvalidPhysicalInputs(testCase)
            testCase.verifyError(@() uav.geometry.worldToArray(diag([1, 1, 2]), [0; 0; -1]), ...
                'uav:geometry:InvalidRotation');
            testCase.verifyError(@() uav.geometry.worldToArray(eye(3), [0; 0; -2]), ...
                'uav:geometry:InvalidDirectionNorm');
            testCase.verifyError(@() uav.geometry.directionCoordinates([1.01; 0; 0]), ...
                'uav:geometry:InvalidDirectionNorm');
            testCase.verifyError(@() uav.array.steeringUPA([1.01; 0], 2, 2, [0.5, 0.5]), ...
                'uav:array:OutsideVisibleDisk');

            [xi, azElRad] = uav.geometry.directionCoordinates([0; 0; 1]);
            a = uav.array.steeringUPA(xi, 2, 2, [0.5, 0.5]);
            testCase.verifyEqual(azElRad, [0; pi / 2], 'AbsTol', 1e-12);
            testCase.verifyEqual(a, ones(4, 1), 'AbsTol', 1e-12);
        end

        % 显式验证矩形UPA使用vx⊗vy顺序，并且每列范数为sqrt(M)。
        function rectangularSteeringUsesVxKronVyOrderAndHasExpectedNorm(testCase)
            a = uav.array.steeringUPA([0.25; -0.5], 2, 3, [0.5, 0.5]);
            expected = [1; -1i; -1; ...
                (1 + 1i) / sqrt(2); (1 - 1i) / sqrt(2); (-1 - 1i) / sqrt(2)];

            testCase.verifyEqual(a, expected, 'AbsTol', 1e-12);
            testCase.verifyEqual(abs(a), ones(6, 1), 'AbsTol', 1e-12);
            testCase.verifyEqual(norm(a), sqrt(6), 'AbsTol', 1e-12);
        end

        % 长轨迹经验协方差应接近配置稳态协方差；零抖动轨迹应严格不变。
        function scenarioHasStationaryAttitudeCovarianceAndZeroJitterBehavior(testCase)
            cfg = config('smoke');
            cfg.array.Nx = 1;
            cfg.array.Ny = 1;
            cfg.sources.interfererAzElDeg = zeros(0, 2);
            cfg.sources.interfererInrDb = [];
            cfg.timeline.K = 20000;
            cfg = uav.config.validate(cfg);

            scenario = uav.scenario.generate(cfg, uav.util.makeStream(13, 1));
            empiricalVariance = var(scenario.attitudeOffsetRad, 0, 2);
            targetVariance = diag(scenario.Qeta);

            testCase.verifyEqual(scenario.Qnu, scenario.Qeta - scenario.F * ...
                scenario.Qeta * scenario.F', 'AbsTol', 1e-14);
            testCase.verifyEqual(empiricalVariance, targetVariance, 'RelTol', 0.35);

            cfg.attitude.stdDeg = [0, 0, 0];
            cfg = uav.config.validate(cfg);
            zeroScenario = uav.scenario.generate(cfg, uav.util.makeStream(13, 1));
            testCase.verifyEqual(zeroScenario.attitudeOffsetRad, zeros(3, cfg.timeline.K), ...
                'AbsTol', 0);
            testCase.verifyEqual(zeroScenario.rotations(:, :, 1), cfg.derived.nominalRotation, ...
                'AbsTol', 1e-12);
            testCase.verifyEqual(zeroScenario.steeringTruth(:, :, 1), cfg.derived.nominalSteering, ...
                'AbsTol', 1e-12);
        end

        % 同一随机子流逐元素复现快拍，并检查truth与算法观测字段分离。
        function periodIsDeterministicAndSeparatesTruthFromAlgorithmInputs(testCase)
            cfg = config('smoke');
            scenarioA = uav.scenario.generate(cfg, uav.util.makeStream(22, 1));
            scenarioB = uav.scenario.generate(cfg, uav.util.makeStream(22, 1));
            obsA = uav.scenario.generatePeriod(cfg, scenarioA, 3, uav.util.makeStream(22, 2));
            obsB = uav.scenario.generatePeriod(cfg, scenarioB, 3, uav.util.makeStream(22, 2));

            testCase.verifyEqual(obsA.X, obsB.X, 'AbsTol', 0);
            testCase.verifyEqual(obsA.Rx, obsB.Rx, 'AbsTol', 0);
            testCase.verifyEqual(obsA.pilot, obsB.pilot, 'AbsTol', 0);
            testCase.verifySize(obsA.X, [cfg.derived.M, cfg.timeline.L]);
            testCase.verifySize(obsA.Rx, [cfg.derived.M, cfg.derived.M]);
            testCase.verifyEqual(obsA.Rx, obsA.X * obsA.X' / cfg.timeline.L, 'AbsTol', 1e-12);
            testCase.verifyFalse(isfield(obsA, 'trueASVs'));
            testCase.verifyFalse(isfield(obsA, 'trueIncm'));
            testCase.verifyTrue(isfield(obsA, 'truth'));
            testCase.verifySize(obsA.truth.trueASVs, [cfg.derived.M, cfg.derived.J + 1]);
            testCase.verifySize(obsA.truth.trueIncm, [cfg.derived.M, cfg.derived.M]);
            testCase.verifyTrue(all(isfinite(obsA.X), 'all'));
            testCase.verifyTrue(all(isfinite(obsA.Rx), 'all'));
        end

        % 验证L<M且J=0时仍能生成有限快拍和正确纯噪声INCM。
        function periodSupportsFewerSnapshotsThanSensorsAndZeroInterferers(testCase)
            cfg = config('smoke');
            cfg.sources.interfererAzElDeg = zeros(0, 2);
            cfg.sources.interfererInrDb = [];
            cfg.timeline.L = 3;
            cfg = uav.config.validate(cfg);
            scenario = uav.scenario.generate(cfg, uav.util.makeStream(29, 1));
            obs = uav.scenario.generatePeriod(cfg, scenario, 1, uav.util.makeStream(29, 2));

            testCase.verifyLessThan(cfg.timeline.L, cfg.derived.M);
            testCase.verifySize(obs.X, [cfg.derived.M, 3]);
            testCase.verifySize(obs.truth.trueASVs, [cfg.derived.M, 1]);
            testCase.verifyEqual(obs.truth.trueIncm, cfg.derived.noisePower * eye(cfg.derived.M), ...
                'AbsTol', 1e-12);
        end

        % 压力测试J>M，确保干扰数不被阵元数或固定示例数量写死。
        function periodSupportsMoreInterferersThanSensors(testCase)
            cfg = config('smoke');
            cfg.array.Nx = 1;
            cfg.array.Ny = 2;
            cfg.sources.interfererAzElDeg = [40, -25; -55, -20; 70, -30];
            cfg.sources.interfererInrDb = [20, 15, 12];
            cfg = uav.config.validate(cfg);
            scenario = uav.scenario.generate(cfg, uav.util.makeStream(31, 1));
            obs = uav.scenario.generatePeriod(cfg, scenario, 1, uav.util.makeStream(31, 2));

            testCase.verifyGreaterThanOrEqual(cfg.derived.J, cfg.derived.M);
            testCase.verifySize(obs.truth.trueASVs, [2, 4]);
            testCase.verifyTrue(all(isfinite(obs.truth.trueIncm), 'all'));
        end

        % 大L下样本协方差应向由源功率和真实流形得到的理论协方差收敛。
        function largeSnapshotScmConvergesToTheoreticalCovariance(testCase)
            cfg = config('smoke');
            cfg.timeline.L = 50000;
            cfg = uav.config.validate(cfg);
            scenario = uav.scenario.generate(cfg, uav.util.makeStream(37, 1));
            obs = uav.scenario.generatePeriod(cfg, scenario, 1, uav.util.makeStream(37, 2));
            A = obs.truth.trueASVs;
            expectedRx = A * diag(obs.truth.sourcePowers) * A' + ...
                cfg.derived.noisePower * eye(cfg.derived.M);

            testCase.verifyEqual(obs.Rx, expectedRx, 'RelTol', 0.08);
        end

        % 检查QPSK单位模、独立性和复高斯噪声的经验功率统计。
        function qpskSignalsAndDerivedNoiseHaveExpectedStatistics(testCase)
            cfg = config('smoke');
            cfg.timeline.L = 50000;
            cfg = uav.config.validate(cfg);
            scenario = uav.scenario.generate(cfg, uav.util.makeStream(43, 1));
            obs = uav.scenario.generatePeriod(cfg, scenario, 1, uav.util.makeStream(43, 2));
            A = obs.truth.trueASVs;
            desired = sqrt(obs.truth.desiredPower) * A(:, 1) * obs.pilot;
            interference = A(:, 2:end) * (sqrt(obs.truth.interfererPowers(:)) .* ...
                obs.truth.interfererSymbols);
            noise = obs.X - desired - interference;
            pilotInterfererCorrelation = mean(obs.pilot .* conj( ...
                obs.truth.interfererSymbols(1, :)));

            testCase.verifyEqual(abs(obs.pilot), ones(1, cfg.timeline.L), 'AbsTol', 1e-12);
            testCase.verifyEqual(abs(obs.truth.interfererSymbols), ...
                ones(cfg.derived.J, cfg.timeline.L), 'AbsTol', 1e-12);
            testCase.verifyLessThan(abs(pilotInterfererCorrelation), 0.02);
            testCase.verifyEqual(mean(abs(noise).^2, 'all'), cfg.derived.noisePower, ...
                'RelTol', 0.03);
            testCase.verifyLessThan(abs(mean(noise.^2, 'all')), 0.01);
        end
    end
end
