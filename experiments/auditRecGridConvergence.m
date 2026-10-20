function differences=auditRecGridConvergence(mainDirectory,denseDirectory)
%AUDITRECGRIDCONVERGENCE 合并同几何的事前网格与事后加密，保留全部方向点。
load(fullfile(mainDirectory,'analysis.mat'),'analysis');a=analysis.summary;
load(fullfile(denseDirectory,'analysis.mat'),'analysis');b=analysis.summary;
load(fullfile(denseDirectory,'manifest.mat'),'manifest');cfg=manifest.cfg;
T=[a(a.block=="numerical"&a.arrayType=="UPA",:);b];
grids=sort(unique(T.gridSize));methods=cfg.experiments.recComprehensive.methodNames;
geometries=unique(T.geometryId);rows=struct([]);
for mi=1:numel(methods)
    for gi=1:numel(geometries)
        for si=1:numel(grids)-1
            low=T(T.method==methods(mi)&T.geometryId==geometries(gi)&T.gridSize==grids(si),:);
            high=T(T.method==methods(mi)&T.geometryId==geometries(gi)&T.gridSize==grids(si+1),:);
            assert(height(low)==1&&height(high)==1,'网格对照不完整或存在重复。');
            row=struct('method',methods(mi),'geometryId',geometries(gi), ...
                'nominalPriorValid',low.nominalPriorValid,'oldGrid',grids(si),'newGrid',grids(si+1), ...
                'oldSinrDb',low.medianSinrDb,'newSinrDb',high.medianSinrDb, ...
                'differenceDb',high.medianSinrDb-low.medianSinrDb, ...
                'absoluteDifferenceDb',abs(high.medianSinrDb-low.medianSinrDb));
            if isempty(rows),rows=row;else,rows(end+1)=row;end %#ok<AGROW>
        end
    end
end
differences=struct2table(rows);writetable(differences,fullfile(denseDirectory,'grid_convergence_differences.csv'));
summary=groupsummary(differences,{'method','oldGrid','newGrid'},{'mean','max'},'absoluteDifferenceDb');
writetable(summary,fullfile(denseDirectory,'grid_convergence_summary.csv'));disp(summary);
f=figure('Visible','off','Position',[100,100,1600,800]);tiledlayout(1,2);
for transition=[401,801]
    nexttile;matrix=nan(numel(methods),numel(geometries));
    for mi=1:numel(methods)
        for gi=1:numel(geometries)
            q=differences(differences.method==methods(mi)&differences.geometryId==geometries(gi)&differences.oldGrid==transition,:);
            matrix(mi,gi)=q.differenceDb;
        end
    end
    imagesc(matrix);colorbar;yticks(1:numel(methods));
    yticklabels({'原始Rec','旧积分static','历史CovFit','网格CovFit'});
    xlabel('原预设几何序号');title(sprintf('%d → %d：SINR变化 (dB)',transition,2*transition-1));
end
set(findall(f,'-property','FontSize'),'FontSize',cfg.plot.fontSize);
savefig(f,fullfile(denseDirectory,'grid_convergence.fig'));
exportgraphics(f,fullfile(denseDirectory,'grid_convergence.png'),'Resolution',cfg.plot.resolutionDpi);
exportgraphics(f,fullfile(denseDirectory,'grid_convergence.pdf'),'ContentType','vector');close(f);
auditSource=fileread(mfilename('fullpath')+".m");
save(fullfile(denseDirectory,'grid_convergence_audit.mat'),'differences','summary','auditSource');
end
