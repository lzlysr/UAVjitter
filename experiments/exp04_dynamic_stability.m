function result = exp04_dynamic_stability(cfg)
%EXP04_DYNAMIC_STABILITY 动态稳定性实验的公开入口。
%   返回逐试次逐周期原始指标、统计摘要和完整结果文件路径。
result = uav.experiment.runExp04(cfg);
end
