function result = exp03_snapshot_sweep(cfg)
%EXP03_SNAPSHOT_SWEEP 有限快拍数扫描的公开入口。
%   入口只转发配置，不在实验脚本中散落快拍数或算法参数。
result = uav.experiment.runExp03(cfg);
end
