function [fig, statistics] = jitterMismatchBoxplot(jitterPoints, rawMismatch, cfg)
%JITTERMISMATCHBOXPLOT 绘制实验2的归一化期望导向失配箱线图。
%   rawMismatch维度为[独立试次×相关周期×抖动扫描点]。每个试次先跨
%   周期取中位数，再用这些试次中位数画一个箱子，不能把周期当独立样本。
%   箱体上下边为Q1/Q3，中线为中位数；须不是置信区间，也不是5%/95%
%   分位线，而是距箱体不超过配置倍数IQR的最远实际样本。所有离群点
%   均画出，不剔除数据。只用原生图元，无需依赖boxplot工具箱接口。

jitterPoints = reshape(jitterPoints, 1, []);
if ~isnumeric(rawMismatch) || ~isreal(rawMismatch) || isempty(rawMismatch) || ...
        ndims(rawMismatch) > 3 || any(~isfinite(rawMismatch), 'all') || ...
        any(rawMismatch < 0, 'all') || ...
        ~isnumeric(jitterPoints) || ~isreal(jitterPoints) || ...
        isempty(jitterPoints) || any(~isfinite(jitterPoints)) || ...
        numel(unique(jitterPoints)) ~= numel(jitterPoints) || ...
        size(rawMismatch, 3) ~= numel(jitterPoints)
    error('uav:plot:InvalidMismatchSamples', ...
        '失配必须为[试次×周期×扫描点]的有限非负实数数组，扫描点须唯一。');
end
style = cfg.plot.exp02Box;
% reshape保留单试次、单扫描点等退化维度，不能用squeeze猜测维度。
samples = reshape(median(rawMismatch, 2), size(rawMismatch, 1), []);
pointCount = numel(jitterPoints);
quartiles = quantile(samples, [0.25, 0.5, 0.75], 1);
iqrValues = quartiles(3, :) - quartiles(1, :);
whiskers = zeros(2, pointCount);
outliers = cell(1, pointCount);
for index = 1:pointCount
    values = samples(:, index);
    lowerFence = quartiles(1, index) - style.whiskerIqrMultiplier * iqrValues(index);
    upperFence = quartiles(3, index) + style.whiskerIqrMultiplier * iqrValues(index);
    inside = values >= lowerFence & values <= upperFence;
    whiskers(:, index) = [min(values(inside)); max(values(inside))];
    outliers{index} = values(~inside);
end
statistics = struct('aggregationUnit', "trial-median-over-periods", ...
    'jitterRmsDeg', jitterPoints, 'trialMedian', samples, ...
    'quartiles', quartiles, 'whiskers', whiskers, 'outliers', {outliers}, ...
    'whiskerIqrMultiplier', style.whiskerIqrMultiplier, ...
    'trialCount', size(samples, 1), 'periodCount', size(rawMismatch, 2));

% 保留定量横轴的实际间距；箱宽只取最小扫描间隔的一定比例。
if pointCount > 1
    spacing = min(diff(sort(jitterPoints)));
else
    spacing = style.singlePointSpacing;
end
halfWidth = spacing * style.boxWidthRatio / 2;
capHalfWidth = halfWidth * style.capWidthRatio;
figureCfg = cfg;
figureCfg.plot.figurePosition = style.figurePosition;
fig = uav.plot.newFigure(figureCfg);
ax = axes(fig);
hold(ax, 'on');
for index = 1:pointCount
    x = jitterPoints(index);
    color = style.colors(mod(index - 1, size(style.colors, 1)) + 1, :);
    % 须与端帽先画，随后用白底空心箱体盖住箱体内部的须线。
    line(ax, [x, x], whiskers(:, index), ...
        'Color', style.whiskerColor, 'LineWidth', style.whiskerLineWidth);
    for endpoint = 1:2
        line(ax, x + [-capHalfWidth, capHalfWidth], ...
            repmat(whiskers(endpoint, index), 1, 2), ...
            'Color', style.whiskerColor, 'LineWidth', style.whiskerLineWidth);
    end
    boxHandle = patch(ax, x + [-halfWidth, halfWidth, halfWidth, -halfWidth], ...
        quartiles([1, 1, 3, 3], index), 'w', ...
        'EdgeColor', color, 'LineWidth', style.boxLineWidth, 'Tag', 'MismatchBox');
    line(ax, x + [-halfWidth, halfWidth], repmat(quartiles(2, index), 1, 2), ...
        'Color', color, 'LineWidth', style.medianLineWidth, ...
        'Tag', 'MismatchMedian', 'Clipping', 'off');
    if index == 1
        % 图例仅说明箱体，引用第一个箱体以避免重复条目。
        legendBox = boxHandle;
    end
    if ~isempty(outliers{index})
        line(ax, repmat(x, size(outliers{index})), outliers{index}, ...
            'LineStyle', 'none', 'Marker', '.', 'Color', color, ...
            'MarkerSize', style.outlierMarkerSize, 'Tag', 'MismatchOutlier');
    end
end
xlim(ax, [min(jitterPoints), max(jitterPoints)] + ...
    [-1, 1] * spacing * style.xPaddingSteps);
if isempty(style.yLimits)
    % 相对范数误差可超过1；自动上限必须覆盖离群点，不能按“归一化”截断。
    limits = ylim(ax);
    ylim(ax, [0, limits(2)]);
else
    ylim(ax, style.yLimits);
end
xlabel(ax, '姿态抖动RMS（°）', 'FontSize', cfg.plot.fontSize);
ylabel(ax, '归一化期望导向失配程度', 'FontSize', cfg.plot.fontSize);
set(ax, 'FontSize', cfg.plot.fontSize, 'LineWidth', style.axesLineWidth, ...
    'Box', 'on', 'TickDir', 'in', 'XGrid', 'off', 'YGrid', 'off');
% 仅保留箱体图例，中位线、须和离群点仍正常绘制。
legend(ax, legendBox, {'箱体：25%–75%分位'}, ...
    'FontSize', cfg.plot.fontSize, 'Location', 'northoutside', ...
    'AutoUpdate', 'off', 'Interpreter', 'none');
% FIG中一并保存绘图样本与统计定义，便于检查箱体和须的来源。
fig.UserData = struct('boxplot', statistics, 'plotStyle', style, ...
    'fontSize', cfg.plot.fontSize, ...
    'description', '箱体：25%–75%分位；中线：中位数；须：IQR规则内实际极值；点：离群样本。');
end
