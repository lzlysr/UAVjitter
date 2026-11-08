function halfWidthDeg = desiredExclusionHalfWidth(cfg)
%DESIREDEXCLUSIONHALFWIDTH 返回当前阵列模式实际使用的期望排除半宽（度）。
%   新配置的ULA使用专用参数；缺少该字段的历史MAT沿用旧公共宽度，
%   以免重新加载旧配置时静默改变实验条件。UPA始终使用原公共字段。
halfWidthDeg = cfg.spectrum.desiredExclusionHalfWidthDeg;
if cfg.array.Nx>1 && cfg.array.Ny==1 && ...
        isfield(cfg.spectrum,'ulaDesiredExclusionHalfWidthDeg')
    halfWidthDeg = cfg.spectrum.ulaDesiredExclusionHalfWidthDeg;
end
end
