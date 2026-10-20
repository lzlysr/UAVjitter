classdef TestUlaGeometryStudy < matlab.unittest.TestCase
    %TESTULAGEOMETRYSTUDY 检验一维移植的流形、真值隔离及几何反例。
    methods(Test)
        function linearManifoldMatchesDegenerateUpa(testCase)
            cfg=config('smoke'); g=cfg.experiments.recUlaGeometry;
            theta=deg2rad([0,-40,55]);
            first=uav.array.steeringULA(theta,g.sensorCount,g.spacingOverLambda);
            second=uav.array.steeringUPA([zeros(size(theta));sin(theta)], ...
                1,g.sensorCount,[g.spacingOverLambda,g.spacingOverLambda]);
            testCase.verifyEqual(first,second,'AbsTol',1e-14);
        end
        function refinedCovfitUsesOnlyAllowedInputs(testCase)
            [cfg,obs,A,Ri]=fixture([0,-40,55],8);
            [Rh,d]=uav.incm.ulaCovfit(obs,0,cfg,25,true);
            obs.truth=struct('directions',NaN,'covariance',NaN);
            cfg.sources.interfererAzElDeg=[];
            cfg.derived.nominalSteering=NaN;
            repeated=uav.incm.ulaCovfit(obs,0,cfg,25,true);
            testCase.verifyEqual(Rh,repeated,'AbsTol',0);
            testCase.verifyFalse(d.usedTruth);
            testCase.verifyLessThan(norm(Rh-Rh','fro'),1e-10);
            testCase.verifyGreaterThan(min(real(eig(Rh))),0);
            [w,~]=uav.beamforming.solve('mvdr',Rh,A(:,1),cfg);
            value=10*log10(cfg.derived.desiredPower*abs(w'*A(:,1))^2/real(w'*Ri*w));
            testCase.verifyGreaterThan(value,19);
            testCase.verifyEqual(w'*A(:,1),1,'AbsTol',1e-10);
        end
        function endfireRecLossExistsWithExactInput(testCase)
            [cfg,obs,A,Ri]=fixture([80,-50,-20],64);
            g=cfg.experiments.recUlaGeometry; count=max(g.rawGridSizes);
            step=pi/count; theta=-pi/2+step/2+step*(0:count-1);
            mask=abs(rad2deg(theta)-80)>25;
            G=uav.array.steeringULA(theta(mask),g.sensorCount,g.spacingOverLambda);
            Rh=uav.incm.guLeshem2012Rec(obs.Rx,G,step*ones(1,sum(mask)));
            [w,d]=uav.beamforming.solve('mvdr',Rh,A(:,1),cfg);
            value=cfg.derived.desiredPower*abs(w'*A(:,1))^2/real(w'*Ri*w);
            upper=cfg.derived.desiredPower*real(A(:,1)'*(Ri\A(:,1)));
            testCase.verifyGreaterThan(10*log10(upper/value),7);
            testCase.verifyEqual(d.loading,0);
        end
        function nearInterferersRemainFiniteWithoutTruthCount(testCase)
            [cfg,obs]=fixture([5,-31,-29],8);
            for refine=[false,true]
                [Rh,d]=uav.incm.ulaCovfit(obs,deg2rad(5),cfg,25,refine);
                testCase.verifyTrue(all(isfinite(Rh),'all'));
                testCase.verifyGreaterThan(min(real(eig(Rh))),0);
                testCase.verifyLessThanOrEqual(numel(d.nodesRad)-d.desiredCount, ...
                    cfg.experiments.recUlaGeometry.maximumCandidateCount);
            end
        end
    end
end

function [cfg,obs,A,Ri]=fixture(angles,L)
%FIXTURE 人工总体输入独立于实验包装，X只提供快拍数分支夹具。
cfg=config('smoke'); g=cfg.experiments.recUlaGeometry;
A=uav.array.steeringULA(deg2rad(angles),g.sensorCount,g.spacingOverLambda);
Ri=cfg.derived.noisePower*eye(g.sensorCount)+A(:,2:end)*diag(cfg.derived.interfererPowers)*A(:,2:end)';
obs=struct('X',zeros(g.sensorCount,L),'Rx',Ri+cfg.derived.desiredPower*A(:,1)*A(:,1)');
end
