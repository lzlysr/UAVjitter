function cfg = applyClassicParameterSelection(cfg, selection)
%APPLYCLASSICPARAMETERSELECTION 将已锁定的经典方法参数应用到新配置。
%   该函数用于把smoke阶段选出的参数原样迁移到preliminary配置。它只
%   修改与四种待筛选方法有关的字段，不改变场景、随机种子、抖动扫描点
%   或蒙特卡洛规模，因而不会发生依据preliminary结果再次调参的情况。

required = {'snapshotCount', 'dlDiagonalLoading', ...
    'espSignalSubspaceDimension', 'wcoDesignJitterStdDeg', ...
    'staticSpectrumDiagonalLoading', ...
    'staticDesiredExclusionHalfWidthDeg', 'staticIncmShrinkage'};
if ~isstruct(selection) || ~all(isfield(selection, required))
    error('uav:experiment:InvalidClassicParameterSelection', ...
        '锁定参数必须包含快拍数以及DL、ESP、WCO和静态INCM的全部字段。');
end

cfg.timeline.L = selection.snapshotCount;
cfg.beamforming.diagonalLoading = selection.dlDiagonalLoading;
cfg.beamforming.classic.esp.signalSubspaceDimension = ...
    selection.espSignalSubspaceDimension;
cfg.beamforming.calibration.designJitterStdDeg = ...
    selection.wcoDesignJitterStdDeg * ones(1, 3);
cfg.spectrum.diagonalLoading = selection.staticSpectrumDiagonalLoading;
cfg.spectrum.desiredExclusionHalfWidthDeg = ...
    selection.staticDesiredExclusionHalfWidthDeg;
cfg.incm.shrinkage = selection.staticIncmShrinkage;
cfg = uav.config.validate(cfg);
end
