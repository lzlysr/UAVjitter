function result=diagnose_upa_integral_grid(profile)
%DIAGNOSE_UPA_INTEGRAL_GRID 独立执行原始Rec的UPA网格收敛诊断，不加入run_all。
%   startup; diagnose_upa_integral_grid("preliminary")
%   数值设置在config.m的experiments.upaGridConvergence中。
if nargin<1,profile="preliminary";end
cfg=config(profile);
result=uav.experiment.runUpaIntegralGridConvergence(cfg);
end
