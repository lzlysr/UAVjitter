function fig = newFigure(cfg)
%NEWFIGURE 按统一配置创建实验图窗。
%   批处理默认Visible='off'，但仍生成可编辑FIG、矢量PDF和PNG。

if cfg.plot.visible
    visibility = 'on';
else
    visibility = 'off';
end
fig = figure('Visible', visibility, ...
    'Position', cfg.plot.figurePosition, 'Color', 'w');
end
