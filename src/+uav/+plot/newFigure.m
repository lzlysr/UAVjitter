function fig = newFigure(cfg)
%NEWFIGURE Create a configured noninteractive experiment figure.

if cfg.plot.visible
    visibility = 'on';
else
    visibility = 'off';
end
fig = figure('Visible', visibility, ...
    'Position', cfg.plot.figurePosition, 'Color', 'w');
end
