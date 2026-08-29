function result = exp06_classic_method_jitter_sweep(cfg)
%EXP06_CLASSIC_METHOD_JITTER_SWEEP 经典方法姿态抖动退化实验公开入口。
%   方法组合、抖动扫描点、试次和周期数全部从cfg读取；该入口不会被
%   run_all自动调用，适合在形成“已完成的初步工作”图表时单独运行。
result = uav.experiment.runClassicJitterSweep(cfg);
end
