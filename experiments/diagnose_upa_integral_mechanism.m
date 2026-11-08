function result=diagnose_upa_integral_mechanism(profile)
%DIAGNOSE_UPA_INTEGRAL_MECHANISM 研究分支的二维积分原理控制实验。
%   startup; result=diagnose_upa_integral_mechanism("preliminary");
%   所有数值设置在config.m中；不改正式实验6/7/8，也不研究processed方法。
if nargin<1,profile="preliminary";end
result=uav.experiment.runUpaIntegralMechanism(config(profile));
end
