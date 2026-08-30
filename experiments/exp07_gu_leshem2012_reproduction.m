function result = exp07_gu_leshem2012_reproduction(cfg)
%EXP07_GU_LESHEM2012_REPRODUCTION 复现原论文例1的Rec比较曲线。
%   本实验独立于UAV主场景，不由run_all自动运行；输出对应论文图1(a)
%   的SNR扫描和图1(c)的快拍数扫描。
result = uav.experiment.runGuLeshem2012Reproduction(cfg);
end
