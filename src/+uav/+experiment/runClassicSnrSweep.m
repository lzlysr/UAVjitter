function result = runClassicSnrSweep(cfg)
%RUNCLASSICSNRSWEEP 在固定快拍数、固定姿态RMS下扫描期望输入SNR。
%   扫描点、方法组合和固定条件均来自config.m的实验8专用字段。
%   SNR相对单阵元噪声功率定义；干扰INR与噪声功率不随SNR改变。
%   与实验6/7共用方法调度、配对观测、统计及文件导出，不加入run_all。
result = uav.experiment.runClassicSweep(cfg, "snr");
end
