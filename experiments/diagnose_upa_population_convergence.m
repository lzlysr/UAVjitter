function result=diagnose_upa_population_convergence(cfg)
%DIAGNOSE_UPA_POPULATION_CONVERGENCE 原始Rec的零抖动理论协方差网格诊断。
%   参数集中于config.m的experiments.upaPopulationConvergence。
%   所有结果写入独立、唯一的北京时间目录，不覆盖任何既有实验文件。
result=uav.experiment.runUpaPopulationConvergence(cfg);
end
