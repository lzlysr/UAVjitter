function paths = plotLinYueSnapshotStudy(matPath)
%PLOTLINYUESNAPSHOTSTUDY 从已保存的配对原始结果出图，不重跑、不挑选试次。
%   第一图展示5°目标场景；第二图同时展示0°/5°候选和功率消融。
%   Lin/Yue标签均指重构核心移植，不代表完整BCS/DSPR算法复现。
loaded=load(matPath,'result');r=loaded.result;cfg=r.cfg;
paths=strings(2,3);
fig=uav.plot.newFigure(cfg);ax=axes(fig,'Position',[.10,.25,.86,.65]);
ji=find(r.jitterStdDeg==max(r.jitterStdDeg),1);
methods=[1,2,3,4,5,22,23,24];
if any(r.methodNames=="le_refined_covfit"),methods=[methods,25,26];end
draw(ax,methods,ji);
title(ax,sprintf('3×3 UPA，三轴姿态RMS=%g°；固定名义期望导向',r.jitterStdDeg(ji)));
legend(ax,cellstr(r.methodLabels(methods)),'Location','southoutside','NumColumns',3,'FontSize',cfg.plot.fontSize);
firstStatus=r.raw.status(:,:,1,ji,23);
firstFailed=nnz(firstStatus==-1);firstSuccess=nnz(firstStatus==1);
annotation(fig,'textbox',[.08,.01,.87,.08],'String', ...
    {'Lin/Yue：功率重构核心＋共同的观测驱动干扰候选细化，非完整论文算法。', ...
    sprintf('Lin在L=%d时原式SCM不可逆；Yue该点%d/%d次非PSD，仅统计%d个成功样本；其余点失败率见CSV。', ...
    r.scanPoints(1),firstFailed,numel(firstStatus),firstSuccess)}, ...
    'EdgeColor','none','FontSize',cfg.plot.fontSize);
paths(1,:)=publish(fig,"lin_yue_"+cfg.profile+"_comparison");

fig=uav.plot.newFigure(cfg);layout=tiledlayout(fig,1,numel(r.jitterStdDeg));
methods=[5,6,10,19,21,22,15];
if any(r.methodNames=="le_refined_covfit"),methods=[methods,25,26];end
for ji=1:numel(r.jitterStdDeg)
    ax=nexttile(layout);draw(ax,methods,ji);
    title(ax,sprintf('三轴姿态RMS=%g°',r.jitterStdDeg(ji)));
end
legend(ax,cellstr(r.methodLabels(methods)),'Location','southoutside','FontSize',cfg.plot.fontSize);
paths(2,:)=publish(fig,"lin_yue_"+cfg.profile+"_causal_ablation");

    function draw(ax,methods,ji)
        % 每条曲线使用不同标记；Oracle固定黑色虚线。所有坐标显式设置对数轴。
        colors=lines(numel(methods));markers={'none','o','s','^','v','d','p','h','x','+'};
        styles={'--','-','-.',':','-','--','-.',':','-','--'};
        hold(ax,'on');
        for index=1:numel(methods)
            color=colors(index,:);if methods(index)==1,color=[0,0,0];end
            plot(ax,r.scanPoints,r.summary{ji}.medianDb(:,methods(index)), ...
                'Color',color,'LineStyle',styles{index},'Marker',markers{index}, ...
                'LineWidth',cfg.plot.lineWidth,'MarkerSize',7);
        end
        set(ax,'XScale','log','XTick',r.scanPoints,'FontSize',cfg.plot.fontSize);
        grid(ax,'on');xlabel(ax,'每周期快拍数 L','FontSize',cfg.plot.fontSize);
        ylabel(ax,'输出SINR中位数（dB）','FontSize',cfg.plot.fontSize);
    end
    function files=publish(fig,name)
        files=strings(1,3);formats=cfg.plot.formats;
        for index=1:numel(formats)
            files(index)=fullfile(cfg.output.figureDir,name+"."+formats(index));
            if formats(index)=="fig"
                savefig(fig,files(index));
            elseif formats(index)=="pdf"
                exportgraphics(fig,files(index),'ContentType','vector');
            else
                exportgraphics(fig,files(index),'Resolution',cfg.plot.resolutionDpi);
            end
        end
        close(fig);
    end
end
