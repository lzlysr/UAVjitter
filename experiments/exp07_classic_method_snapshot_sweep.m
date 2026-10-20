function result = exp07_classic_method_snapshot_sweep(cfg)
%EXP07_CLASSIC_METHOD_SNAPSHOT_SWEEP 六方法快拍数扫描的公开入口。
%   与实验6共享经典方法调度，不修改原实验3，也不自动加入run_all。
result = uav.experiment.runClassicSnapshotSweep(cfg);
end
