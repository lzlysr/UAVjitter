function [A,D,d] = leSharedCandidates(obs,original,cfg)
%LESHAREDCANDIDATES 将Le全分区细化器接到本次相同局部峰，仅细化干扰。
%   这是机制对照，不是Le完整方法：保留期望字典、候选列数和顺序，
%   不进行细化后去重，便于隔离局部有界搜索与全分区搜索的区别。
%   最终期望参考由调用方固定；本函数不读取真值或真实干扰数。
[~,D,nodes]=uav.spectrum.localPeakCandidates(obs.Rx,original,cfg);
Rc=(obs.Rx+obs.Rx')/2+cfg.incm.staticCovfit.loadingRatio* ...
    real(trace(obs.Rx))/cfg.derived.M*eye(cfg.derived.M);
[nodes(D+1:end,:),d]=uav.spectrum.refineCaponNodesLe(nodes(D+1:end,:),Rc,false,cfg);
A=uav.array.steeringUPA(nodes.',cfg.array.Nx,cfg.array.Ny,cfg.array.spacingOverLambda);
d.desiredCount=D;d.usedTruth=false;
end
