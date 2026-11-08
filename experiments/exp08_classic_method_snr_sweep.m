function result = exp08_classic_method_snr_sweep(cfg)
%EXP08_CLASSIC_METHOD_SNR_SWEEP 经典方法输入SNR扫描的公开入口。
%   默认扫描-10:5:30 dB、32快拍、三轴各5度RMS抖动；参数统一读config.m。
result = uav.experiment.runClassicSnrSweep(cfg);
end
