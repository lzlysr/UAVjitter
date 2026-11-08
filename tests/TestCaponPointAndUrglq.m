classdef TestCaponPointAndUrglq < matlab.unittest.TestCase
    % 核验公式、求积域、真值隔离以及实验6/7/8的小规模接入，不跑正式规模。
    methods(TestMethodSetup)
        function paths(t)
            root=fileparts(fileparts(mfilename('fullpath')));
            for path={root,fullfile(root,'src'),fullfile(root,'experiments')}
                t.applyFixture(matlab.unittest.fixtures.PathFixture(path{1}));
            end
        end
    end
    methods(Test)
        function configurationAndInvalidControls(t)
            c=config("smoke");names=["incm_capon_point","incm_urglq"];
            for field={'classicJitterMethodNames','classicSnapshotMethodNames','classicSnrMethodNames'}
                t.verifyTrue(all(ismember(names,c.experiments.(field{1}))));
            end
            invalid=c;invalid.incm.urglq.quadratureOrder=0;
            t.verifyError(@()uav.config.validate(invalid),'uav:config:InvalidAlgorithmParameter');
            invalid=c;invalid.incm.urglq.upaAzElHalfWidthDeg=[8,-1];
            t.verifyError(@()uav.config.validate(invalid),'uav:config:InvalidAlgorithmParameter');
            invalid=c;invalid.incm.caponPoint.powerLoadingRatio=-1;
            t.verifyError(@()uav.config.validate(invalid),'uav:config:InvalidAlgorithmParameter');
            invalid=c;invalid.incm.caponCandidates.peakThresholdRatio=2;
            t.verifyError(@()uav.config.validate(invalid),'uav:config:InvalidAlgorithmParameter');
            invalid=rmfield(c.incm,'urglq');missing=c;missing.incm=invalid;
            t.verifyError(@()uav.config.validate(missing),'uav:config:InvalidExperimentParameter');
        end
        function pointFormulaAndSourceTruthIsolation(t)
            c=smallConfig([3,3]);obs=observation(c);safe=uav.experiment.algorithmObservation(obs);
            [R,~,d]=uav.incm.reconstruct("incm_capon_point",safe,[],c);
            A=d.candidateManifold;p=1./real(sum(conj(A).*(safe.Rx\A),1));
            expected=(A.*p)*A'+d.noisePower*eye(c.derived.M);
            t.verifyEqual(R,expected,'AbsTol',1e-9);
            t.verifyEqual(d.caponPowers,p.','AbsTol',1e-10);
            t.verifyFalse(d.usedTruth);t.verifyFalse(d.powerFittingApplied);
            t.verifyEqual(d.powerSpectrum.loading,0);
            % 保持相同观测，删除配置中的真实干扰列表；算法不能借此改变输出。
            altered=c;altered.sources.interfererAzElDeg=zeros(0,2);altered.sources.interfererInrDb=[];
            altered=uav.config.validate(altered);
            other=uav.incm.reconstruct("incm_capon_point",safe,[],altered);
            t.verifyEqual(R,other);
            other=uav.incm.reconstruct("incm_urglq",safe,[],altered);
            original=uav.incm.reconstruct("incm_urglq",safe,[],c);
            t.verifyEqual(other,original);
            t.verifyError(@()uav.incm.reconstruct("incm_capon_point",obs,[],c),'uav:incm:InvalidObservation');
        end
        function ulaThreePointRuleAndOverlappingIntervals(t)
            c=smallConfig([9,1]);center=deg2rad(40);half=deg2rad(c.incm.urglq.ulaHalfWidthDeg);
            [nodes,w,d]=uav.spectrum.urglqQuadrature([sin(center),0],c);
            expectedAngle=center+half*[-sqrt(3/5);0;sqrt(3/5)];
            t.verifyEqual(asin(nodes(:,1)),expectedAngle,'AbsTol',1e-14);
            t.verifyEqual(w,half*[5/9;8/9;5/9],'AbsTol',1e-14);
            t.verifyEqual(d.measure,"dtheta");
            [~,w,d]=uav.spectrum.urglqQuadrature([sind(40),0;sind(45),0],c);
            t.verifyEqual(sum(w),deg2rad(21),'AbsTol',1e-14);
            t.verifyEqual(d.domainCellCount,1);
        end
        function upaSolidAngleUnionAndPeriodicBoundary(t)
            c=smallConfig([3,3]);c.spectrum.desiredExclusionHalfWidthDeg=0;c=uav.config.validate(c);
            h=deg2rad(c.incm.urglq.upaAzElHalfWidthDeg);el=deg2rad(-25);
            candidates=[cos(el)*cosd(40),cos(el)*sind(40)];
            [nodes,w,d]=uav.spectrum.urglqQuadrature(candidates,c);
            area=2*h(1)*(sin(el+h(2))-sin(el-h(2)));
            t.verifySize(nodes,[9,2]);t.verifyEqual(sum(w),area,'AbsTol',1e-10);
            t.verifyEqual(d.measure,"solid_angle");
            next=[cos(el)*cosd(45),cos(el)*sind(45)];
            [~,w]=uav.spectrum.urglqQuadrature([candidates;next;candidates],c);
            t.verifyEqual(sum(w),deg2rad(21)*(sin(el+h(2))-sin(el-h(2))),'AbsTol',1e-10);
            nearHorizon=deg2rad(-1);seam=[cos(nearHorizon)*cosd(179),cos(nearHorizon)*sind(179)];
            [nodes,w]=uav.spectrum.urglqQuadrature(seam,c);
            t.verifyTrue(all(abs(atan2d(nodes(:,2),nodes(:,1)))>170));
            t.verifyEqual(sum(w),2*h(1)*(sin(0)-sin(nearHorizon-h(2))),'AbsTol',1e-10);
            [nodes,w]=uav.spectrum.urglqQuadrature(zeros(0,2),c);
            t.verifySize(nodes,[0,2]);t.verifyEmpty(w);
        end
        function urglqIsWeightedIntegralWithoutProjection(t)
            for shape={[3,3],[9,1]}
                c=smallConfig(shape{1});obs=observation(c);safe=uav.experiment.algorithmObservation(obs);
                [R,~,d]=uav.incm.reconstruct("incm_urglq",safe,[],c);
                A=uav.array.steeringUPA(d.nodes.',c.array.Nx,c.array.Ny,c.array.spacingOverLambda);
                p=1./real(sum(conj(A).*(safe.Rx\A),1));
                expected=(A.*(d.weights.'.*p))*A';
                t.verifyNotEmpty(d.weights);t.verifyEqual(R,expected,'AbsTol',1e-9);
                t.verifyFalse(d.desiredProjectionApplied);t.verifyFalse(d.desiredSteeringCorrectionApplied);
                t.verifyFalse(d.completePaperReproduction);t.verifyFalse(d.usedTruth);
                t.verifyEqual(R,R','AbsTol',1e-10);
                t.verifyGreaterThanOrEqual(min(real(eig(R))),-1e-9*norm(R,2));
                record=uav.experiment.runClassicMethods(obs,[],c,["oracle","incm_urglq"]);
                expectedWeight=uav.beamforming.solve("mvdr",R,c.derived.nominalSteeringReference,c);
                t.verifyEqual(record.weights(:,2),expectedWeight,'AbsTol',1e-10);
                t.verifyLessThanOrEqual(record.sinrLinear(2),record.sinrLinear(1)*(1+1e-9));
            end
        end
        function singularScmIsExplicitlyNotApplicable(t)
            c=smallConfig([3,3]);c.timeline.L=8;c=uav.config.validate(c);obs=observation(c);
            names=["incm_capon_point","incm_urglq"];
            record=uav.experiment.runClassicMethods(obs,[],c,names);
            t.verifyEqual(record.methodStatus,["not-applicable","not-applicable"]);
            t.verifyTrue(all(isnan(record.sinrLinear)));
            c.incm.caponPoint.powerLoadingRatio=1e-3;c.incm.urglq.powerLoadingRatio=1e-3;
            record=uav.experiment.runClassicMethods(obs,[],c,names);
            t.verifyEqual(record.methodStatus,["ok","ok"]);
            t.verifyTrue(all(isfinite(record.sinrLinear)));
            for k=1:2,t.verifyGreaterThan(record.diagnostics{k}.reconstruction.powerSpectrum.loading,0);end
        end
        function bothArraysAndThreeSweepsExport(t)
            folder=t.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            for shape={[3,3],[9,1]}
                c=smallConfig(shape{1});names=["oracle","incm_capon_point","incm_urglq"];
                c.experiments.classicJitterMethodNames=names;
                c.experiments.classicSnapshotMethodNames=names;c.experiments.classicSnrMethodNames=names;
                c.experiments.jitterSweepDeg=0;c.experiments.snapshotSweep=32;c.experiments.snrSweepDb=10;
                c.experiments.classicSnrJitterStdDeg=[0,0,0];
                c.monteCarlo.nTrials=1;c.experiments.bootstrapReplicates=10;
                for field={'jitterSweep','snapshotSweep','snrSweep'},c.experiments.periodCounts.(field{1})=1;end
                c.plot.visible=false;c.plot.closeAfterExport=true;
                for field={'rawDir','figureDir','tableDir','logDir'}
                    c.output.(field{1})=fullfile(folder.Folder,field{1});
                end
                for scan=["jitter","snapshots","snr"]
                    r=uav.experiment.runClassicSweep(c,scan);
                    t.verifyEqual(r.status,"complete");t.verifyEqual(r.methodNames,names);
                    t.verifyTrue(all(isfinite(r.raw.sinrLinear),'all'));
                    t.verifyEqual(r.raw.snapshotSignatures(:,:,:,1,:),r.raw.snapshotSignatures(:,:,:,2,:));
                    t.verifyEqual(r.raw.snapshotSignatures(:,:,:,1,:),r.raw.snapshotSignatures(:,:,:,3,:));
                    for field={'mat','fig','pdf','png','csv','log'},t.verifyTrue(isfile(r.artifacts.(field{1})));end
                end
            end
        end
    end
end

function c=smallConfig(shape)
% 测试夹具仅缩小工作量；物理条件与算法默认值保持不变。
c=config("smoke");c.array.Nx=shape(1);c.array.Ny=shape(2);
c.spectrum.ulaGridSize=201;c.attitude.stdDeg=[0,0,0];
c=uav.config.validate(c);
end
function obs=observation(c)
s=uav.scenario.generate(c,uav.util.makeStream(c.random.baseSeed,1));
obs=uav.scenario.generatePeriod(c,s,1,uav.util.makeStream(c.random.baseSeed,2));
end
