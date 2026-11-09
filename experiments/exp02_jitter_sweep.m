function result = exp02_jitter_sweep(cfg)
%EXP02_JITTER_SWEEP 姿态抖动强度扫描的公开入口。
%   所有扫描点与蒙特卡洛规模均从cfg读取，正式实现位于runExp02。
%   失配图直接使用箱线图，沿用原实验结果文件名，不另建重绘入口或后缀。
result = uav.experiment.runExp02(cfg);
end
