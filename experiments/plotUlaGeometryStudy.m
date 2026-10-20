function paths=plotUlaGeometryStudy(profile,keepOpen)
%PLOTULAGEOMETRYSTUDY 导出一维几何诊断曲线，保留全部有利和不利点。
arguments
    profile (1,1) string = "preliminary"
    keepOpen (1,1) logical = true
end
cfg=config(profile); g=cfg.experiments.recUlaGeometry;
root=fileparts(fileparts(mfilename('fullpath')));
directory=fullfile(root,g.outputDirectory,profile);
loaded=load(fullfile(directory,'ula_geometry.mat'),'results'); r=loaded.results; t=r.records;
directory=fullfile(directory,'figures'); if ~isfolder(directory), mkdir(directory); end
families=unique(string({r.geometries.family}),'stable');
names=["原始Rec（25°）","原始Rec（7°）","网格CovFit（25°）", ...
    "网格CovFit（7°）","细化CovFit（25°）","细化CovFit（7°）"];
methods=["raw_rec","raw_rec","covfit_grid","covfit_grid","covfit_refined","covfit_refined"];
halves=repmat(g.halfWidthsDeg,1,3); colors=lines(numel(names));
styles=["-","--",":","-.","-","--"]; markers=["o","s","^","d","v","p"];
xLabels=["期望ULA角度（°）","全部源一维角度平移（°）", ...
    "第一干扰相对期望的角度偏移（°）","两个干扰的角度间隔（°）","已测UPA世界几何编号"];
paths=strings(0,1);
for familyIndex=1:numel(families)
    indices=find(string({r.geometries.family})==families(familyIndex));
    x=[r.geometries(indices).parameter];
    for L=[Inf,g.snapshotCounts]
        fig=uav.plot.newFigure(cfg); ax=axes(fig); hold(ax,'on');
        upper=arrayfun(@(index) t.oracleSinrDb(find(t.geometryIndex==index,1)),indices);
        plot(ax,x,upper,'k--','DisplayName','Oracle','LineWidth',cfg.plot.lineWidth);
        for methodIndex=1:numel(names)
            if methods(methodIndex)=="raw_rec" && L<g.sensorCount, continue; end
            values=nan(size(indices));
            for k=1:numel(indices)
                select=t.geometryIndex==indices(k) & t.L==L & t.method==methods(methodIndex) & ...
                    t.halfWidthDeg==halves(methodIndex);
                if isinf(L) && methods(methodIndex)=="raw_rec"
                    select=select & t.rawGridSize==max(g.rawGridSizes);
                end
                values(k)=10*log10(median(t.sinrLinear(select)));
            end
            plot(ax,x,values,'Color',colors(methodIndex,:),'LineStyle',styles(methodIndex), ...
                'Marker',markers(methodIndex),'DisplayName',names(methodIndex), ...
                'LineWidth',cfg.plot.lineWidth,'MarkerSize',cfg.plot.markerSize);
        end
        xlabel(ax,xLabels(familyIndex),'FontSize',cfg.plot.fontSize);
        ylabel(ax,'输出SINR（dB）','FontSize',cfg.plot.fontSize);
        if isinf(L), description="精确协方差"; else, description="L="+L+"，"+g.trialCount+"试次中位数"; end
        title(ax,"9阵元ULA："+description,'FontSize',cfg.plot.fontSize);
        set(ax,'FontSize',cfg.plot.fontSize); grid(ax,'on'); box(ax,'on');
        legend(ax,'show','Location','bestoutside','FontSize',cfg.plot.fontSize);
        stem=fullfile(directory,families(familyIndex)+"_L"+L);
        savefig(fig,stem+".fig");
        exportgraphics(fig,stem+".png",'Resolution',cfg.plot.resolutionDpi);
        exportgraphics(fig,stem+".pdf",'ContentType','vector');
        paths=[paths;stem+[".fig";".png";".pdf"]]; %#ok<AGROW>
        if ~keepOpen, close(fig); end
    end
end
end
