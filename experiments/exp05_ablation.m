function result = exp05_ablation(cfg)
%EXP05_ABLATION 时序稀疏INCM模块消融的公开入口。
%   变体顺序和单开关约束由config与runVariants共同固定。
result = uav.experiment.runExp05(cfg);
end
