classdef TestUlaAzimuthAndIntegrationGrid < matlab.unittest.TestCase
    %TESTULAAZIMUTHANDINTEGRATIONGRID 核查ULA方位角语义及四方法共用稠密网格。
    methods(TestMethodSetup)
        function addPaths(testCase)
            root=fileparts(fileparts(mfilename('fullpath')));
            for path={root,fullfile(root,'src')}
                testCase.applyFixture(matlab.unittest.fixtures.PathFixture(path{1}));
            end
        end
    end
    methods(Test)
        function nominalAzimuthUsesFirstColumnAndHeightIsIgnored(testCase)
            c=ulaConfig();angles=[c.sources.desiredAzElDeg(1),c.sources.interfererAzElDeg(:,1).'];
            expected=exp(1i*pi*(0:8)'*sind(angles));
            testCase.verifyEqual(c.derived.nominalSteering,expected,'AbsTol',1e-12);
            testCase.verifyEqual(rad2deg(c.derived.ulaGeometry.nominalAnglesRad),angles,'AbsTol',1e-12);
            altered=c;altered.sources.desiredAzElDeg(2)=75;
            altered.sources.interfererAzElDeg(:,2)=[0;120];
            altered=uav.config.validate(altered);
            testCase.verifyEqual(c.derived.worldUnitVectors,altered.derived.worldUnitVectors);
            [scene,obs]=observation(c);[otherScene,otherObs]=observation(altered);
            testCase.verifyEqual(scene.steeringTruth,otherScene.steeringTruth);
            testCase.verifyEqual(obs.X,otherObs.X);
            testCase.verifyTrue(c.derived.ulaGeometry.heightIgnored);
        end
        function fullRotationStillDrivesActualUlaResponse(testCase)
            c=ulaConfig();scene=observation(c);
            for k=1:c.timeline.K
                direction=scene.rotations(:,:,k).'*c.derived.worldUnitVectors;
                expected=exp(1i*pi*(0:8)'*direction(1,:));
                testCase.verifyEqual(scene.steeringTruth(:,:,k),expected,'AbsTol',1e-12);
                testCase.verifyEqual(scene.ulaAzimuthRad(:,k),asin(direction(1,:)).','AbsTol',1e-12);
            end
            testCase.verifyNotEqual(scene.steeringTruth(:,:,1),scene.steeringTruth(:,:,2));
        end
        function denseIntegralIsIndependentFromCandidateGrid(testCase)
            c=ulaConfig();g=c.derived.ulaIntegrationGrid;
            testCase.verifyEqual(size(c.derived.spectrumGrid.manifold,2),c.spectrum.gridSize);
            testCase.verifyEqual(size(g.manifold,2),2001);
            testCase.verifyEqual(sum(g.weights),pi,'AbsTol',1e-12);
            Rx=diag(1:9);keep=uav.spectrum.admissibleDirections(g.directionCosines.',c);
            A=g.manifold(:,keep);p=1./real(sum(conj(A).*(Rx\A),1));
            expected=(A.*(reshape(g.weights(keep),1,[]).*p))*A';
            [actual,~,d]=uav.incm.reconstruct("incm_integral_origin",struct('Rx',Rx),[],c);
            testCase.verifyEqual(actual,expected,'AbsTol',1e-11);
            testCase.verifyEqual(d.integrationGridSize,2001);
            [~,~,d]=uav.incm.reconstruct("incm_integral_processed",struct('Rx',Rx),[],c);
            testCase.verifyEqual(d.integrationGridSize,2001);
            testCase.verifyEqual(numel(d.powerMap),2001);
            % 处理后积分应等同于只把候选网格设为同密度的旧求积实现。
            dense=c;dense.spectrum.gridSize=2001;
            dense.spectrum=rmfield(dense.spectrum,'ulaGridSize');dense=uav.config.validate(dense);
            testCase.verifyEqual(uav.incm.reconstruct("incm_integral_processed",struct('Rx',Rx),[],c), ...
                uav.incm.reconstruct("incm_integral_processed",struct('Rx',Rx),[],dense),'AbsTol',1e-12);
        end
        function changingComparisonGridDoesNotChangeTemporal(testCase)
            c=ulaConfig();[~,obs]=observation(c);safe=uav.experiment.algorithmObservation(obs);
            changed=c;changed.spectrum.ulaGridSize=801;changed=uav.config.validate(changed);
            testCase.verifyEqual(c.derived.spectrumGrid,changed.derived.spectrumGrid);
            R=uav.incm.reconstruct("temporal",safe,[],c);
            testCase.verifyEqual(R,uav.incm.reconstruct("temporal",safe,[],changed),'AbsTol',0);
        end
        function bothCovfitsUseTheSameDenseGridThroughAllEntrypoints(testCase)
            c=ulaConfig();[~,obs]=observation(c);safe=uav.experiment.algorithmObservation(obs);
            for G=[801,2001]
                cfg=c;cfg.spectrum.ulaGridSize=G;cfg=uav.config.validate(cfg);
                % 等价参照：把基础谱本身设为G点，并删除独立网格参数。
                reference=cfg;reference.spectrum.gridSize=G;
                reference.spectrum=rmfield(reference.spectrum,'ulaGridSize');
                reference=uav.config.validate(reference);
                for method=["incm_covfit","incm_covfit_refined"]
                    if method=="incm_covfit",entry=@uav.incm.upaCovfit;else,entry=@uav.incm.localCovfit;end
                    [actual,d]=entry(safe,cfg);
                    expected=entry(safe,reference);
                    testCase.verifyEqual(actual,expected,'AbsTol',1e-12);
                    testCase.verifyEqual(actual,uav.incm.reconstruct(method,safe,[],cfg),'AbsTol',1e-12);
                    testCase.verifyEqual(d.candidateGridSize,G);
                    testCase.verifyEqual(d.candidateGridStep,pi/G,'AbsTol',1e-15);
                    if method=="incm_covfit_refined"
                        testCase.verifyEqual(d.refinement.radiusRad,cfg.incm.localCovfit.radiusBins*pi/G,'AbsTol',1e-15);
                    end
                end
            end
            % 选网格仅作用于局部副本，用户配置的基础谱缓存不被覆盖。
            testCase.verifyEqual(size(c.derived.spectrumGrid.manifold,2),c.spectrum.gridSize);
        end
        function classicDispatcherUsesOneGridForAllFourMethods(testCase)
            c=ulaConfig();[~,obs]=observation(c);
            methods=["incm_integral_origin","incm_integral_processed","incm_covfit","incm_covfit_refined"];
            r=uav.experiment.runClassicMethods(obs,[],c,methods);
            for k=1:4
                d=r.diagnostics{k}.reconstruction;
                if k<=2,actual=d.integrationGridSize;else,actual=d.candidateGridSize;end
                testCase.verifyEqual(actual,c.spectrum.ulaGridSize);
                testCase.verifyEqual(r.snapshotSignatures(k,:),r.snapshotSignatures(1,:));
                testCase.verifyEqual(r.diagnostics{k}.beamforming.method,"mvdr");
            end
        end
        function denseQuadratureConvergesForDefaultZeroJitterGeometry(testCase)
            c=ulaConfig();c.attitude.stdDeg(:)=0;
            % 仅检查当前几何/INR的数值收敛，不把该阈值当作所有场景的保证。
            for snr=[-10,10,30]
                c.signal.snrDb=snr;c=uav.config.validate(c);[~,obs]=observation(c);
                a=obs.truth.trueASVs(:,1);
                safe=struct('Rx',obs.truth.trueIncm+c.derived.desiredPower*(a*a'));
                fine=c;fine.spectrum.ulaGridSize=4001;fine=uav.config.validate(fine);
                for method=["incm_integral_origin","incm_integral_processed"]
                    first=uav.incm.reconstruct(method,safe,[],c);
                    second=uav.incm.reconstruct(method,safe,[],fine);
                    w1=uav.beamforming.solve("mvdr",first,a,c);
                    w2=uav.beamforming.solve("mvdr",second,a,fine);
                    m1=uav.metrics.evaluate(w1,obs.truth,c);m2=uav.metrics.evaluate(w2,obs.truth,fine);
                    testCase.verifyLessThan(abs(m1.sinrDb-m2.sinrDb),0.1);
                end
            end
        end
        function upaKeepsThreeDimensionalModelAndOriginalIntegral(testCase)
            c=config("smoke");c.array.Nx=3;c.array.Ny=3;c=uav.config.validate(c);
            testCase.verifyFalse(isfield(c.derived,'ulaIntegrationGrid'));
            testCase.verifyEqual(c.derived.worldUnitVectors(:,1),[cosd(35);0;-sind(35)],'AbsTol',1e-12);
            changed=c;changed.spectrum.ulaGridSize=801;
            changed.array.ulaDirectionConvention="world_az_el_projection";changed=uav.config.validate(changed);
            testCase.verifyEqual(c.derived,changed.derived);
            for method=["incm_integral_origin","incm_integral_processed"]
                testCase.verifyEqual(uav.incm.reconstruct(method,struct('Rx',diag(1:9)),[],c), ...
                    uav.incm.reconstruct(method,struct('Rx',diag(1:9)),[],changed),'AbsTol',0);
            end
            c.sources.desiredAzElDeg(2)=0;
            testCase.verifyError(@()uav.config.validate(c),'uav:config:SourceNotLowerHemisphere');
        end
        function historicalConfigurationKeepsProjectionAndCoarseGrid(testCase)
            c=ulaConfig();c.array=rmfield(c.array,'ulaDirectionConvention');
            c.spectrum=rmfield(c.spectrum,'ulaGridSize');c=uav.config.validate(c);
            testCase.verifyEqual(rad2deg(c.derived.ulaGeometry.nominalAnglesRad(1)),55,'AbsTol',1e-12);
            testCase.verifyFalse(isfield(c.derived,'ulaIntegrationGrid'));
            testCase.verifyEqual(c,uav.config.validate(c));
        end
        function historicalGridFieldMigratesWithoutChangingResults(testCase)
            c=ulaConfig();old=c;
            old.spectrum=rmfield(old.spectrum,'ulaGridSize');
            old.spectrum.ulaIntegrationGridSize=c.spectrum.ulaGridSize;
            testCase.verifyEqual(uav.config.validate(old),c);
            % 两个字段同时存在时以用户显式设置的新字段为准。
            both=c;both.spectrum.ulaIntegrationGridSize=801;
            testCase.verifyEqual(uav.config.validate(both),c);
        end
        function invalidNewSettingsAreRejected(testCase)
            c=ulaConfig();bad=c;bad.sources.desiredAzElDeg(1)=91;
            testCase.verifyError(@()uav.config.validate(bad),'uav:config:InvalidDirection');
            for G=[0,2,4,NaN]
                bad=c;bad.spectrum.ulaGridSize=G;
                testCase.verifyError(@()uav.config.validate(bad),'uav:config:InvalidPositiveInteger');
            end
        end
    end
end

function c=ulaConfig()
%ULACONFIG 显式采用x轴ULA，避免依赖用户当前默认阵列形状。
c=config("smoke");c.array.Nx=9;c.array.Ny=1;c=uav.config.validate(c);
end

function [scene,obs]=observation(c)
%OBSERVATION 相同子流用于高度角不变性及观测配对检查。
scene=uav.scenario.generate(c,uav.util.makeStream(c.random.baseSeed,81));
obs=uav.scenario.generatePeriod(c,scene,1,uav.util.makeStream(c.random.baseSeed,82));
end
