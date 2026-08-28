function result = exp01_model_validation(cfg)
%EXP01_MODEL_VALIDATION 建模与机理验证的公开入口。
%   cfg应由config(profile)生成。入口不保存任何数值，只把配置交给正式
%   实现；返回结构包含原始结果、统计摘要和六种格式的产物路径。
result = uav.experiment.runExp01(cfg);
end
