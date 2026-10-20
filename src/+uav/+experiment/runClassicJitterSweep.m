function result = runClassicJitterSweep(cfg)
%RUNCLASSICJITTERSWEEP 保留经典方法姿态抖动扫描的原公开接口。
%   方法组合由classicJitterMethodNames控制，每周期快拍数固定为timeline.L。
result = uav.experiment.runClassicSweep(cfg, "jitter");
end
