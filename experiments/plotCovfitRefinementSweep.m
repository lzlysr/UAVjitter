function fig = plotCovfitRefinementSweep(result)
%PLOTCOVFITREFINEMENTSWEEP 绘制四个单因素轴，不把搜索半径错误施加给Le全域方法。
cfg=result.cfg;fig=uav.plot.newFigure(cfg);layout=tiledlayout(fig,2,2);
colors=lines(numel(result.methodNames));markers={'none','o','s','^','v','d','p'};
styles={'--','-','-.',':','--','-.',':'};
for ai=1:numel(result.axisNames)
    ax=nexttile(layout);hold(ax,'on');
    for mi=1:numel(result.methodNames)
        color=colors(mi,:);if result.methodNames(mi)=="oracle",color=[0,0,0];end
        plot(ax,result.scanPoints{ai},result.summary{ai}.medianDb(:,mi), ...
            'Color',color,'LineStyle',styles{mi},'Marker',markers{mi}, ...
            'LineWidth',cfg.plot.lineWidth,'MarkerSize',cfg.plot.markerSize);
    end
    set(ax,'XScale','log','XTick',sort(result.scanPoints{ai}),'FontSize',cfg.plot.fontSize);
    grid(ax,'on');xlabel(ax,result.axisLabels(ai),'FontSize',cfg.plot.fontSize);
    ylabel(ax,'输出SINR中位数（dB）','FontSize',cfg.plot.fontSize);
end
lg=legend(ax,cellstr(result.methodLabels),'FontSize',cfg.plot.fontSize,'NumColumns',3);
lg.Layout.Tile='south';
title(layout,sprintf('3×3 UPA，L=%d，三轴RMS=%g°；固定名义期望导向', ...
    cfg.experiments.linYueStudy.refinementSweep.snapshotCount, ...
    cfg.experiments.linYueStudy.refinementSweep.jitterStdDeg),'FontSize',cfg.plot.fontSize);
end
