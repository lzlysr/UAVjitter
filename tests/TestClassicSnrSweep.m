classdef TestClassicSnrSweep < matlab.unittest.TestCase
    %TESTCLASSICSNRSWEEP 检查实验8的固定条件、配对观测、统计和输出。
    methods (TestMethodSetup)
        function addProjectPaths(testCase)
            root=fileparts(fileparts(mfilename('fullpath')));
            for path={root,fullfile(root,'src'),fullfile(root,'experiments')}
                testCase.applyFixture(matlab.unittest.fixtures.PathFixture(path{1}));
            end
        end
    end
    methods (Test)
        function defaultsAndValidation(testCase)
            cfg=config("smoke");
            testCase.verifyEqual(cfg.experiments.snrSweepDb,-10:5:30);
            testCase.verifyEqual(cfg.experiments.classicSnrSnapshotCount,32);
            testCase.verifyEqual(cfg.experiments.classicSnrJitterStdDeg,[5,5,5]);
            cfg.experiments.classicSnrMethodNames=["oracle","CovFit"];
            checked=uav.config.validate(cfg);
            testCase.verifyEqual(checked.experiments.classicSnrMethodNames,["oracle","incm_covfit"]);
            cfg.experiments.classicSnrMethodNames=["oracle","oracle"];
            testCase.verifyError(@() uav.config.validate(cfg),'uav:config:InvalidExperimentParameter');
            cfg.experiments.classicSnrMethodNames="unknown";
            testCase.verifyError(@() uav.config.validate(cfg),'uav:config:InvalidExperimentParameter');
        end
        function rejectsInvalidPhysicalAndScanParameters(testCase)
            original=config("smoke");
            cfg=original;cfg.experiments.snrSweepDb=[0,NaN];
            testCase.verifyError(@() uav.config.validate(cfg),'uav:config:InvalidVector');
            cfg=original;cfg.experiments.classicSnrSnapshotCount=1.5;
            testCase.verifyError(@() uav.config.validate(cfg),'uav:config:InvalidPositiveInteger');
            cfg=original;cfg.experiments.classicSnrJitterStdDeg=[0,-1,0];
            testCase.verifyError(@() uav.config.validate(cfg),'uav:config:InvalidExperimentParameter');
        end
        function historicalConfigsRemainValidWithoutExperiment8(testCase)
            cfg=config("smoke");
            cfg.experiments=rmfield(cfg.experiments,{'snrSweepDb','classicSnrSnapshotCount', ...
                'classicSnrJitterStdDeg','classicSnrMethodNames'});
            cfg.experiments.periodCounts=rmfield(cfg.experiments.periodCounts,'snrSweep');
            checked=uav.config.validate(cfg);
            testCase.verifyFalse(isfield(checked.experiments,'snrSweepDb'));
        end
        function bothArraysKeepFixedConditionsAndExport(testCase)
            folder=testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            for shape={[3,3],[9,1]}
                cfg=config("smoke");cfg.array.Nx=shape{1}(1);cfg.array.Ny=shape{1}(2);
                % 故意令全局条件与实验8不同，以验证专用配置实际覆盖生效。
                cfg.timeline.L=8;cfg.attitude.stdDeg=[2,3,4];
                cfg.monteCarlo.nTrials=1;cfg.experiments.periodCounts.snrSweep=1;
                cfg.experiments.snrSweepDb=[-10,30];
                cfg.experiments.bootstrapReplicates=10;
                cfg.plot.visible=false;cfg.plot.closeAfterExport=true;
                cfg.output.rawDir=fullfile(folder.Folder,'raw');
                cfg.output.figureDir=fullfile(folder.Folder,'figures');
                cfg.output.tableDir=fullfile(folder.Folder,'tables');
                cfg.output.logDir=fullfile(folder.Folder,'logs');
                r=exp08_classic_method_snr_sweep(cfg);
                testCase.verifyEqual(r.status,"complete");
                testCase.verifyEqual(r.cfg.timeline.L,32);
                testCase.verifyEqual(r.cfg.attitude.stdDeg,[5,5,5]);
                testCase.verifyEqual(cfg.timeline.L,8);
                testCase.verifyEqual(cfg.attitude.stdDeg,[2,3,4]);
                testCase.verifyEqual(r.cfg.sources.interfererInrDb,cfg.sources.interfererInrDb);
                testCase.verifyEqual(r.cfg.signal.noisePower,cfg.signal.noisePower);
                testCase.verifyEqual(r.methodNames,cfg.experiments.classicSnrMethodNames);
                testCase.verifyTrue(all(isfinite(r.raw.sinrLinear),'all'));
                % 非零抖动应实际传入模型；不同SNR点不强制使用同一姿态实现。
                testCase.verifyTrue(all(r.raw.desiredSteeringMismatchNormalized>0,'all'));
                for m=2:numel(r.methodNames)
                    testCase.verifyEqual(r.raw.snapshotSignatures(:,:,:,1,:),r.raw.snapshotSignatures(:,:,:,m,:));
                end
                for name={'mat','fig','pdf','png','csv','log'}
                    testCase.verifyTrue(isfile(r.artifacts.(name{1})));
                end
                data=readtable(r.artifacts.csv);
                testCase.verifyTrue(ismember('inputSnrDb',data.Properties.VariableNames));
                saved=load(r.artifacts.mat,'result');
                testCase.verifyEqual(saved.result.scanPoints,[-10,30]);
                testCase.verifyNotEqual(r.seedScheme{1}.scenarioSubstreams, ...
                    uav.experiment.substream(r.cfg,7,1,1,0,"scenario"));
            end
        end
        function zeroJitterRemainsAnExplicitControl(testCase)
            cfg=config("smoke");
            folder=testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            for name={'rawDir','figureDir','tableDir','logDir'}
                cfg.output.(name{1})=fullfile(folder.Folder,name{1});
            end
            % 人工零抖动夹具保留理论斜率检查，不修改正式5度默认值。
            cfg.experiments.classicSnrJitterStdDeg=[0,0,0];
            cfg.experiments.classicSnrMethodNames="oracle";
            cfg.experiments.snrSweepDb=[-10,30];
            cfg.monteCarlo.nTrials=1;cfg.experiments.periodCounts.snrSweep=1;
            cfg.experiments.bootstrapReplicates=10;
            cfg.plot.visible=false;cfg.plot.closeAfterExport=true;
            r=exp08_classic_method_snr_sweep(cfg);
            testCase.verifyEqual(r.cfg.attitude.stdDeg,[0,0,0]);
            testCase.verifyEqual(diff(r.summary.sinr.medianDb),40,'AbsTol',1e-8);
        end
    end
end
