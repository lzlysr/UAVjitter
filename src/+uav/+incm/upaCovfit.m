function [RinHat, diagnostics] = upaCovfit(obs, cfg)
%UPACOVFIT 共用细化版前端和功率拟合，只关闭连续干扰方向细化。
%   参数统一读取incm.localCovfit；最终仍使用固定名义期望导向。
[RinHat,diagnostics]=uav.incm.localCovfit(obs,cfg,false);
end
