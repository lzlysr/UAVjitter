classdef TestAutomaticUla < matlab.unittest.TestCase
    %TESTAUTOMATICULA 核查x轴自动降维、原始Rec、拟合及经典实验导出。
    % 人工测试夹具可以覆盖网格/规模，但不据此宣称任一方法性能领先。
    methods (TestMethodSetup)
        function addProjectPath(testCase)
            root=fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,'src')));
        end
    end
    methods (Test)
        function gridAndManifoldAreOneDimensional(testCase)
            cfg = ulaConfig(); g=cfg.derived.spectrumGrid;
            testCase.verifyEqual(cfg.derived.arrayType,"ULA");
            testCase.verifyEqual(g.dimension,1);
            testCase.verifyEqual(g.measure,"dtheta");
            testCase.verifySize(g.xiX,[1,cfg.spectrum.gridSize]);
            testCase.verifyEqual(sum(g.weights),pi,'AbsTol',1e-12);
            expected=exp(1i*2*pi*cfg.array.spacingOverLambda(1)* ...
                (0:cfg.array.Nx-1)'*sin(g.thetaRad));
            testCase.verifyEqual(g.manifold,expected,'AbsTol',1e-12);
            testCase.verifyEqual(cfg.derived.ulaGeometry.interferersInsideDesiredSector,[false,false]);
            testCase.verifyEqual(uav.config.validate(cfg),cfg);
        end
        function upaGridAndCaponRemainUnchanged(testCase)
            cfg=config('smoke');cfg.array.Nx=3;cfg.array.Ny=3;
            cfg=uav.config.validate(cfg);g=cfg.derived.spectrumGrid;
            axis=linspace(-1,1,cfg.spectrum.gridSize);[x,y]=meshgrid(axis,axis);
            testCase.verifyEqual(g.xiX,x);testCase.verifyEqual(g.xiY,y);
            testCase.verifyEqual(g.cellArea,(axis(2)-axis(1))^2);
            testCase.verifyEqual(g.manifold,uav.array.steeringUPA( ...
                [x(g.visibleMask)';y(g.visibleMask)'],3,3,cfg.array.spacingOverLambda));
            Rx=diag(1:9);[p,d]=uav.spectrum.loadedCapon(Rx,cfg);
            R=Rx+cfg.spectrum.diagonalLoading*trace(Rx)/9*eye(9);
            expected=1./real(sum(conj(g.manifold).*(R\g.manifold),1));
            testCase.verifyEqual(p(g.visibleLinearIndices)',expected,'AbsTol',1e-12);
            testCase.verifyEqual(d.dimension,2);
        end
        function exclusionDependsOnlyOnObservableProjection(testCase)
            cfg=ulaConfig();s=cfg.derived.nominalDirectionCosines(1,1);
            nodes=[s,0;s,0.2;sin(deg2rad(-40)),0];
            testCase.verifyEqual(uav.spectrum.admissibleDirections(nodes,cfg),[false;false;true]);
        end
        function rawRecMatchesIndependentAngularSum(testCase)
            cfg=ulaConfig();g=cfg.derived.ulaIntegrationGrid;Rx=diag(1:9);
            keep=uav.spectrum.admissibleDirections(g.directionCosines.',cfg);
            A=g.manifold(:,keep);
            expected=zeros(9);
            for q=1:size(A,2)
                a=A(:,q);expected=expected+g.cellArea/(a'*(Rx\a))*(a*a');
            end
            [actual,~,d]=uav.incm.reconstruct("incm_integral_origin",struct('Rx',Rx),[],cfg);
            testCase.verifyEqual(actual,expected,'AbsTol',1e-11);
            testCase.verifyEqual(d.loading,0);testCase.verifyFalse(d.usedPseudoInverse);
            testCase.verifyEqual(d.measure,"dtheta");
        end
        function insufficientSnapshotsAreNotSilentlyLoaded(testCase)
            cfg=ulaConfig();cfg.timeline.L=8;cfg=uav.config.validate(cfg);
            obs=observation(cfg);
            r=uav.experiment.runClassicMethods(obs,[],cfg, ...
                ["oracle","incm_integral_origin","incm_integral_processed","incm_covfit","incm_covfit_refined"]);
            testCase.verifyEqual(r.methodStatus(2),"not-applicable");
            testCase.verifyTrue(isnan(r.sinrLinear(2)));
            testCase.verifyTrue(all(isfinite(r.sinrLinear([1,3:5]))));
            testCase.verifyNotEmpty(r.notApplicableReason(2));
            testCase.verifyEqual(r.snapshotSignatures,repmat(r.snapshotSignatures(1,:),5,1));
        end
        function refinedFitAndTemporalAreFinite(testCase)
            cfg=ulaConfig();cfg.timeline.L=64;cfg.spectrum.gridSize=101;
            % 排除宽度改为人工5度夹具，以便两个干扰能参与算法检查。
            cfg.spectrum.ulaDesiredExclusionHalfWidthDeg=5;
            cfg=uav.config.validate(cfg);obs=observation(cfg);
            r=uav.experiment.runClassicMethods(obs,[],cfg, ...
                ["oracle","incm_integral_origin","incm_integral_processed","incm_covfit","incm_covfit_refined"]);
            testCase.verifyTrue(all(r.methodStatus=="ok"));
            testCase.verifyTrue(all(isfinite(r.sinrLinear)));
            testCase.verifyGreaterThanOrEqual(r.sinrLinear(1)+1e-9,max(r.sinrLinear));
            for m=2:5
                testCase.verifyLessThan(abs(r.weights(:,m)'*cfg.derived.nominalSteeringReference-1),1e-10);
            end
            d=r.diagnostics{5}.reconstruction;
            testCase.verifyEqual(d.interfererCandidateNodes(:,2),zeros(size(d.interfererCandidateNodes,1),1));
            testCase.verifyEqual(numel(unique(d.interfererCandidateNodes(:,1))),size(d.interfererCandidateNodes,1));
            testCase.verifyLessThan(d.relativeKktResidual,cfg.incm.localCovfit.kktTolerance);
            safe=uav.experiment.algorithmObservation(obs);state=[];
            for k=1:2
                [R,state,diag]=uav.incm.reconstruct("temporal",safe,state,cfg);
                testCase.verifyTrue(all(isfinite(R),'all'));
                testCase.verifyGreaterThanOrEqual(min(real(eig(R))),-1e-10);
                testCase.verifyEqual(diag.spectrum.dimension,1);
                testCase.verifyEqual(diag.alignment.shiftBins(1),0);
            end
        end
        function adaptiveAngularWeightsConserveMass(testCase)
            cfg=ulaConfig();g=cfg.derived.spectrumGrid;
            [p,s]=uav.spectrum.loadedCapon(eye(9),cfg);
            mask=reshape(uav.spectrum.admissibleDirections(g.directionCosines.',cfg),size(p));
            peaks=find(mask,1);
            q=uav.spectrum.adaptiveQuadrature(p,s,mask,peaks,cfg);
            testCase.verifyEqual(sum(q.weights),nnz(mask)*g.step,'AbsTol',1e-12);
            testCase.verifyTrue(all(uav.spectrum.admissibleDirections(q.nodes,cfg)));
            testCase.verifyEqual(q.nodes(:,2),zeros(size(q.nodes,1),1));
        end
        function emptyInterferenceAndFixedRadiusStayWellDefined(testCase)
            cfg=ulaConfig();cfg.sources.interfererAzElDeg=zeros(0,2);
            cfg.sources.interfererInrDb=zeros(1,0);cfg.attitude.stdDeg(:)=0;
            for L=[8,64]
                cfg.timeline.L=L;cfg=uav.config.validate(cfg);obs=observation(cfg);
                safe=uav.experiment.algorithmObservation(obs);
                for name=["incm_integral_processed","temporal"]
                    R=uav.incm.reconstruct(name,safe,[],cfg);
                    testCase.verifyTrue(all(isfinite(R),'all'));
                end
                for name=["incm_covfit","incm_covfit_refined"]
                    r=uav.experiment.runClassicMethods(obs,[],cfg,name);
                    testCase.verifyTrue(isfinite(r.sinrLinear));
                end
                e=uav.beamforming.calibrateEpsilon(cfg);
                testCase.verifyEqual(e,sqrt(cfg.beamforming.calibration.ulaSquaredRadiusScale*9));
            end
            cfg.beamforming.calibration.ulaSquaredRadiusScale=1;
            testCase.verifyError(@()uav.config.validate(cfg),'uav:config:InvalidAlgorithmParameter');
        end
        function selectedClassicSweepsExportUlaAndMissingRec(testCase)
            cfg=ulaConfig();cfg.monteCarlo.nTrials=1;
            cfg.experiments.bootstrapReplicates=10;
            cfg.experiments.periodCounts.jitterSweep=1;cfg.experiments.periodCounts.snapshotSweep=1;
            cfg.experiments.jitterSweepDeg=[0,2];cfg.experiments.snapshotSweep=[8,64];
            names=["oracle","dl_mvdr","esp","wco","incm_integral_origin","incm_integral_processed","incm_covfit","incm_covfit_refined"];
            cfg.experiments.classicJitterMethodNames=names;cfg.experiments.classicSnapshotMethodNames=names;
            cfg.plot.visible=false;cfg.plot.closeAfterExport=true;
            root=fullfile(cfg.output.rootDir,'results');
            for pair={"rawDir","raw";"figureDir","figures";"tableDir","tables";"logDir","logs"}'
                cfg.output.(pair{1})=fullfile(root,pair{2},'automatic_ula_smoke');
            end
            a=uav.experiment.runClassicJitterSweep(cfg);
            b=uav.experiment.runClassicSnapshotSweep(cfg);
            testCase.verifyEqual(b.experimentName,"exp07_classic_method_snapshot_sweep_ULA_9x1");
            for result={a,b}
                r=result{1};testCase.verifyEqual(r.status,"complete");
                testCase.verifyTrue(contains(r.experimentName,"ULA_9x1"));
                testCase.verifyTrue(isfile(r.artifacts.mat));testCase.verifyTrue(isfile(r.artifacts.png));
                testCase.verifyEqual(r.methodNames,names);
                testCase.verifyEqual(r.methodLabels(5:8),names(5:8));
                t=readtable(r.artifacts.csv);testCase.verifyTrue(ismember('notApplicableReason',t.Properties.VariableNames));
            end
            testCase.verifyEqual(b.summary.sinr.validCounts(1,5),0);
            testCase.verifyEqual(b.raw.methodStatus(1,1,1,5),"not-applicable");
            testCase.verifyEqual(b.calibration.epsilonS,sqrt(0.3*cfg.derived.M),'AbsTol',1e-12);
        end
    end
end

function cfg=ulaConfig()
%ULACONFIG 仅改变阵列形状；其他默认物理与算法参数沿用main。
cfg=config('smoke');cfg.array.Nx=9;cfg.array.Ny=1;
cfg=uav.config.validate(cfg);
end

function obs=observation(cfg)
%OBSERVATION 固定独立子流，生成同一份观测供所有方法使用。
scene=uav.scenario.generate(cfg,uav.util.makeStream(cfg.random.baseSeed,31));
obs=uav.scenario.generatePeriod(cfg,scene,1,uav.util.makeStream(cfg.random.baseSeed,32));
end
