function result = runClassicSnapshotSweep(cfg)
%RUNCLASSICSNAPSHOTSWEEP 经典方法在固定姿态抖动RMS下扫描每周期快拍数。
%   方法组合由classicSnapshotMethodNames控制，扫描点来自snapshotSweep。
%   三轴RMS沿用attitude.stdDeg，方法在同一扫描点内共享全部观测。
result = uav.experiment.runClassicSweep(cfg, "snapshots");
end
