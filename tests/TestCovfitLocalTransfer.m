classdef TestCovfitLocalTransfer < matlab.unittest.TestCase
    %TESTCOVFITLOCALTRANSFER 仅核查移植、改名及实验6/7的小规模接入，不跑全量实验。
    methods(TestMethodSetup)
        function paths(testCase)
            root=fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,'src')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,'experiments')));
        end
    end
    methods(Test)
        function renamedMethodsAndLegacyAlias(testCase)
            cfg=config('smoke');
            testCase.verifyTrue(all(ismember(["CovFit","CovFit_local"],cfg.experiments.classicJitterMethodNames)));
            cfg.experiments.classicSnapshotMethodNames=["static_incm_mvdr","CovFit_local"];
            cfg=uav.config.validate(cfg);
            testCase.verifyEqual(cfg.experiments.classicSnapshotMethodNames,["CovFit","CovFit_local"]);
            cfg.experiments.classicSnapshotMethodNames=["CovFit","static_incm_mvdr"];
            testCase.verifyError(@()uav.config.validate(cfg),'uav:config:InvalidExperimentParameter');
        end
        function newConfigurationRejectsInvalidRadius(testCase)
            cfg=config('smoke');cfg.incm.localCovfit.radiusBins=0;
            testCase.verifyError(@()uav.config.validate(cfg),'uav:config:InvalidAlgorithmParameter');
        end
        function pairedArchivedObservationsMatchResearch(testCase)
            root=fileparts(fileparts(mfilename('fullpath')));base=config('smoke');
            research='D:/workspace/UAVjitter-worktrees/lin-yue-upa-snapshot-study/results/raw/lin_yue_refinement_2026-10-21';
            snapshot=fullfile(research,'lin_yue_snapshot_study_smoke.mat');
            parameter=fullfile(research,'covfit_refinement_parameter_sweep_smoke.mat');
            testCase.assumeTrue(isfile(snapshot)&&isfile(parameter),'需要本机研究分支的smoke归档作为移植夹具。');
            rows={};
            for fixture=1:3
                if fixture<=2
                    loaded=load(snapshot,'result');r=loaded.result;
                    ip=1;if fixture==2,ip=numel(r.scanPoints);end
                    ji=find(r.jitterStdDeg==5);seed=r.seedScheme{ip,ji};L=r.scanPoints(ip);
                    expected=reshape(r.raw.sinrLinear(1,:,ip,ji,[5,22]),[],2);
                else
                    loaded=load(parameter,'result');r=loaded.result;L=32;seed=r.seedScheme;
                    ip=find(r.scanPoints{1}==2);
                    expected=reshape(r.raw.sinrLinear(1,:,ip,1,[2,3]),[],2);
                end
                cfg=r.cfg;cfg.timeline.L=L;cfg.attitude.stdDeg(:)=5;
                cfg.incm.localCovfit=base.incm.localCovfit;cfg=uav.config.validate(cfg);
                scene=uav.scenario.generate(cfg,uav.util.makeStream(cfg.random.baseSeed,seed.scenarioSubstreams(1)));
                for k=1:2
                    obs=uav.scenario.generatePeriod(cfg,scene,k, ...
                        uav.util.makeStream(cfg.random.baseSeed,seed.snapshotSubstreams(1,k)));
                    record=uav.experiment.runClassicMethods(obs,[],cfg,["CovFit","CovFit_local"]);
                    testCase.verifyEqual(record.sinrLinear,expected(k,:),'RelTol',1e-10);
                    testCase.verifyEqual(record.diagnostics{2}.beamforming.method,"mvdr");
                    testCase.verifyFalse(record.diagnostics{2}.reconstruction.usedTruth);
                    testCase.verifyLessThan(record.diagnostics{2}.reconstruction.relativeKktResidual,1e-8);
                    testCase.verifyLessThan(abs(record.weights(:,2)'*cfg.derived.nominalSteeringReference-1),1e-10);
                    old=uav.experiment.runClassicMethods(obs,[],cfg,"static_incm_mvdr");
                    testCase.verifyEqual(old.weights,record.weights(:,1),'AbsTol',1e-12);
                    for mi=1:2
                        rows{end+1}=table(L,k,record.methodNames(mi),10*log10(expected(k,mi)), ...
                            10*log10(record.sinrLinear(mi)),abs(record.sinrLinear(mi)-expected(k,mi)), ...
                            'VariableNames',{'snapshotCount','period','method','researchSinrDb','mainSinrDb','absoluteLinearError'}); %#ok<AGROW>
                    end
                end
            end
            comparison=vertcat(rows{:});disp(comparison);
            folder=fullfile(root,'results','tables','covfit_local_transfer_smoke');
            if ~isfolder(folder),mkdir(folder);end
            writetable(comparison,fullfile(folder,'research_main_parity.csv'));
        end
        function bothClassicExperimentsExportSelectedMethods(testCase)
            root=fileparts(fileparts(mfilename('fullpath')));cfg=config('smoke');
            % 人工轻量夹具只核查接入；图不是20×20研究曲线的重复统计估计。
            cfg.monteCarlo.nTrials=1;cfg.spectrum.gridSize=51;
            cfg.experiments.periodCounts.jitterSweep=2;cfg.experiments.periodCounts.snapshotSweep=2;
            cfg.experiments.jitterSweepDeg=[0,5];cfg.experiments.snapshotSweep=[8,32];
            cfg.experiments.bootstrapReplicates=10;
            cfg.experiments.classicJitterMethodNames=["oracle","CovFit","CovFit_local"];
            cfg.experiments.classicSnapshotMethodNames=cfg.experiments.classicJitterMethodNames;
            cfg.plot.visible=false;cfg.plot.closeAfterExport=true;
            for item={'raw','figure','table','log'}
                field=[item{1} 'Dir'];folder=item{1};if strcmp(folder,'figure'),folder='figures';end
                if strcmp(folder,'table'),folder='tables';end
                if strcmp(folder,'log'),folder='logs';end
                cfg.output.(field)=fullfile(root,'results',folder,'covfit_local_transfer_smoke');
            end
            a=exp06_classic_method_jitter_sweep(cfg);b=exp07_classic_method_snapshot_sweep(cfg);
            for result={a,b}
                r=result{1};testCase.verifyEqual(r.status,"complete");
                testCase.verifyEqual(r.methodNames,["oracle","CovFit","CovFit_local"]);
                testCase.verifyTrue(all(isfinite(r.raw.sinrLinear),'all'));
                testCase.verifyTrue(isfile(r.artifacts.fig));testCase.verifyTrue(isfile(r.artifacts.mat));
                testCase.verifyEqual(r.raw.snapshotSignatures(:,:,:,1,:),r.raw.snapshotSignatures(:,:,:,3,:));
            end
        end
    end
end
