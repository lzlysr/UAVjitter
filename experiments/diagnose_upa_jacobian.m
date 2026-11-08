function result=diagnose_upa_jacobian(profile)
%DIAGNOSE_UPA_JACOBIAN 研究分支专用，不改变正式实验6/7/8或原始Rec定义。
%   startup; result=diagnose_upa_jacobian("preliminary");
%   所有数值设置位于config.m的experiments.upaJacobianStudy。
if nargin<1,profile="preliminary";end
result=uav.experiment.runUpaJacobianStudy(config(profile));
end
