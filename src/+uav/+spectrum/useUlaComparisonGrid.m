function cfg = useUlaComparisonGrid(cfg)
%USEULACOMPARISONGRID 为ULA四种静态重构选择同一高密度角网格。
%   网格点数统一读取ulaGridSize，内部缓存名仍为ulaIntegrationGrid。
%   只修改本次算法调用的配置副本；UPA、旧MAT及temporal的网格不变。
%   CovFit直接入口和统一重构入口均可调用，重复调用不改变结果。
if isfield(cfg,'derived') && isfield(cfg.derived,'ulaIntegrationGrid')
    cfg.derived.spectrumGrid=cfg.derived.ulaIntegrationGrid;
    cfg.spectrum.gridSize=cfg.spectrum.ulaGridSize;
end
end
