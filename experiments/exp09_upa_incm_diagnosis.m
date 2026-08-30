function result = exp09_upa_incm_diagnosis(cfg)
%EXP09_UPA_INCM_DIAGNOSIS 运行UPA二维INCM诊断与协方差拟合实验。
%   本入口独立于run_all，smoke结果只用于验证计算与导出链路。
result = uav.experiment.runUpaIncmDiagnosis(cfg);
end
