function paths = plotRecGeometryStudy(profile,keepOpen)
%PLOTRECGEOMETRYSTUDY 绘制MATLAB探索性几何诊断图，不对扫描点作有利筛选。
%   总体Rec显示已计算的最大网格，其余方法显示各自固定谱网格的结果。
%   有限样本线为独立试次线性SINR中位数；不将扫描点连线解释为时间轨迹。
arguments
    profile (1,1) string = "smoke"
    keepOpen (1,1) logical = true
end
cfg = config(profile);
g = cfg.experiments.recGeometry;
root = fileparts(fileparts(mfilename('fullpath')));
directory = fullfile(root,g.outputDirectory,profile);
loaded = load(fullfile(directory,'geometry_study.mat'),'results');
r = loaded.results; t = r.records;
directory = fullfile(directory,'figures');
if ~isfolder(directory), mkdir(directory); end
families = unique(string({r.geometries.family}),'stable');
xLabels = ["默认几何到正交几何的插值系数", "共同方位旋转（°）", ...
    "第一干扰相对期望的方位偏移（°）", "期望高度角（°）"];
labels = ["Oracle", "原始Rec（25°）", "原始Rec（7°）", "旧static", ...
    "历史CovFit", "细化CovFit（25°）", "细化CovFit（7°）"];
styles = ["--","-","--",":","-.","-","--"];
markers = ["none","o","s","^","d","v","p"];
colors = [0,0,0;lines(numel(r.methodNames)-1)];
paths = strings(0,1);
for familyIndex = 1:numel(families)
    family = families(familyIndex);
    indices = find(string({r.geometries.family})==family);
    for L = [Inf,g.snapshotCounts]
        fig = uav.plot.newFigure(cfg); ax = axes(fig); hold(ax,'on');
        for methodIndex = 1:numel(r.methodNames)
            method = r.methodNames(methodIndex);
            if startsWith(method,"raw_rec") && L<r.cfg.derived.M, continue; end
            values = nan(size(indices));
            for pointIndex = 1:numel(indices)
                selected = t.geometryIndex==indices(pointIndex) & t.L==L & t.method==method;
                if isinf(L) && startsWith(method,"raw_rec")
                    selected = selected & t.rawGridSize==max(g.rawGridSizes);
                end
                values(pointIndex) = 10*log10(median(t.sinrLinear(selected)));
            end
            plot(ax,[r.geometries(indices).parameter],values, ...
                'LineStyle',styles(methodIndex),'Marker',markers(methodIndex), ...
                'Color',colors(methodIndex,:),'DisplayName',labels(methodIndex), ...
                'LineWidth',cfg.plot.lineWidth,'MarkerSize',cfg.plot.markerSize);
        end
        xlabel(ax,xLabels(familyIndex),'FontSize',cfg.plot.fontSize);
        if isinf(L), description = "精确协方差"; else, description = "L="+L+"，"+g.trialCount+"试次中位数"; end
        title(ax,"零抖动："+description,'FontSize',cfg.plot.fontSize);
        ylabel(ax,'输出SINR（dB）','FontSize',cfg.plot.fontSize);
        set(ax,'FontSize',cfg.plot.fontSize); grid(ax,'on'); box(ax,'on');
        legend(ax,'show','Location','bestoutside','FontSize',cfg.plot.fontSize);
        stem = fullfile(directory,family+"_L"+L);
        savefig(fig,stem+".fig");
        exportgraphics(fig,stem+".png",'Resolution',cfg.plot.resolutionDpi);
        exportgraphics(fig,stem+".pdf",'ContentType','vector');
        paths = [paths;stem+[".fig";".png";".pdf"]]; %#ok<AGROW>
        if ~keepOpen, close(fig); end
    end
end
end
