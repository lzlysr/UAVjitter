function result = exp02_jitter_sweep(cfg)
%EXP02_JITTER_SWEEP 姿态抖动强度扫描的公开入口。
%   所有扫描点与蒙特卡洛规模均从cfg读取，正式实现位于runExp02。
result = uav.experiment.runExp02(cfg);
end
