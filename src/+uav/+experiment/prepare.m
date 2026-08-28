function [cfg, calibration] = prepare(cfg)
%PREPARE 在单个实验开始前校验配置并且只标定一次epsilonS。
%   runExp01~05各调用一次本函数。标定子流固定，不随扫描点、试次或周期
%   变化，从而避免使用当前试次真值重新选择鲁棒半径。

cfg = uav.config.validate(cfg);
stream = uav.util.makeStream(cfg.random.baseSeed, ...
    cfg.beamforming.calibration.substream);
[epsilonS, diagnostics] = uav.beamforming.calibrateEpsilon(cfg, stream);
% 后续所有静态/时序INCM方法从同一个cfg读取相同epsilonS。
cfg.beamforming.epsilonS = epsilonS;
calibration = struct( ...
    'epsilonS', epsilonS, ...
    'diagnostics', diagnostics, ...
    'substream', cfg.beamforming.calibration.substream, ...
    'scope', "fixed-once-per-experiment-configuration");
end
