function plotRecComprehensiveStudy(directory)
%PLOTRECCOMPREHENSIVESTUDY 原生MATLAB诊断图；不以精选几何替代完整分布。
%   这是研究诊断图，并非投稿级统计图排版。主比较的同时区间另存CSV。
load(fullfile(directory,'analysis.mat'),'analysis');load(fullfile(directory,'manifest.mat'),'manifest');
T=analysis.summary; g=manifest.cfg.experiments.recComprehensive;
folder=fullfile(directory,'figures');if ~isfolder(folder),mkdir(folder);end
plotSource=fileread(mfilename('fullpath')+".m");save(fullfile(folder,'plot_source.mat'),'plotSource');
labels=["原始Rec","旧积分static","历史CovFit","网格CovFit"];
styles=["-o","--s","-.^",":d"]; colors=lines(numel(labels));
for array=g.arrayTypes
    base=T(T.arrayType==array & T.block=="snapshots",:);
    f=figure('Visible','off','Position',[100,100,1600,800]);tiledlayout(1,2);
    for valid=[true,false]
        nexttile;hold on;selection=base(base.nominalPriorValid==valid,:);
        if isempty(selection),axis off;text(.1,.5,'该子集无几何样本');continue;end
        for mi=1:numel(g.methodNames)
            values=nan(size(g.snapshotCounts));
            for li=1:numel(values)
                q=selection(selection.method==g.methodNames(mi)&selection.L==g.snapshotCounts(li),:);
                values(li)=mean(q.medianLossDb);
            end
            plot(g.snapshotCounts,values,styles(mi),'Color',colors(mi,:),'DisplayName',labels(mi));
        end
        set(gca,'XScale','log');xlabel('快拍数L');ylabel('平均Oracle损失 (dB)');
        title(sprintf('%s：分区先验有效=%d',array,valid));legend('Location','southoutside','NumColumns',2);grid on;
    end
    finish(f,array+"_snapshots");
    f=figure('Visible','off','Position',[100,100,1600,800]);tiledlayout(1,2);
    for valid=[true,false]
        nexttile;hold on;selection=base(base.nominalPriorValid==valid,:);
        if isempty(selection),axis off;text(.1,.5,'该子集无几何样本');continue;end
        for mi=1:numel(g.methodNames)
            values=arrayfun(@(L)mean(selection.medianSinrDb(selection.method==g.methodNames(mi)&selection.L==L)),g.snapshotCounts);
            plot(g.snapshotCounts,values,styles(mi),'Color',colors(mi,:),'DisplayName',labels(mi));
        end
        oracle=arrayfun(@(L)mean(selection.oracleMedianDb(selection.method==g.methodNames(1)&selection.L==L)),g.snapshotCounts);
        plot(g.snapshotCounts,oracle,'k--','DisplayName','理想上界');
        set(gca,'XScale','log');xlabel('快拍数L');ylabel('几何中位数均值 (dB)');
        title(sprintf('%s：SINR，先验有效=%d',array,valid));legend('Location','southoutside','NumColumns',2);grid on;
    end
    finish(f,array+"_sinr");
    population=T(T.arrayType==array&T.block=="population",:);
    ids=unique(population.geometryId); matrix=nan(numel(labels),numel(ids));
    for mi=1:numel(labels)
        for gi=1:numel(ids)
            q=population(population.method==g.methodNames(mi)&population.geometryId==ids(gi),:);
            matrix(mi,gi)=q.medianLossDb;
        end
    end
    f=figure('Visible','off','Position',[100,100,1400,450]);imagesc(matrix);colorbar;
    yticks(1:numel(labels));yticklabels(labels);xlabel('完整几何清单顺序（先随机，后规则网格及锚点）');
    ylabel('方法');title(array+"：总体协方差输入的Oracle损失 (dB)，不裁掉极端值");
    finish(f,array+"_population_all");
    f=figure('Visible','off','Position',[100,100,1200,800]);tiledlayout(2,2);
    for mi=1:numel(labels)
        nexttile;q=population(population.method==g.methodNames(mi),:);
        scatter(q.maxDesiredCoherence,q.medianLossDb,15,double(q.nominalPriorValid),'filled');
        xlabel('最大期望—干扰流形相干性');ylabel('Oracle损失 (dB)');title(labels(mi));grid on;
    end
    finish(f,array+"_coherence");
    f=figure('Visible','off','Position',[100,100,1200,800]);tiledlayout(2,2);
    for mi=1:numel(labels)
        nexttile;q=T(T.arrayType==array&T.block=="power"&T.method==g.methodNames(mi)&T.L==32,:);
        matrix=nan(numel(g.powerSnrDb),numel(g.powerInrDb));
        for si=1:size(matrix,1)
            for ii=1:size(matrix,2)
                matrix(si,ii)=mean(q.medianLossDb(q.snrDb==g.powerSnrDb(si)&q.inrDb==g.powerInrDb(ii)));
            end
        end
        imagesc(g.powerInrDb,g.powerSnrDb,matrix);axis xy;colorbar;
        xlabel('单干扰INR (dB)');ylabel('期望SNR (dB)');title([labels(mi);"L=32，平均损失 (dB)"]);
    end
    finish(f,array+"_power");
    f=figure('Visible','off','Position',[100,100,1600,800]);tiledlayout(1,2);
    for block=["sector","jitter"]
        nexttile;hold on;
        if block=="sector",field='halfWidthDeg';x=g.halfWidthsDeg;xtext='期望区排除半宽 (度)';
        else,field='jitterStdDeg';x=g.jitterStdDeg;xtext='各轴姿态RMS (度)';end
        q=T(T.arrayType==array&T.block==block&T.L==32,:);
        for mi=1:numel(labels)
            y=arrayfun(@(v)mean(q.medianLossDb(q.(field)==v&q.method==g.methodNames(mi))),x);
            plot(x,y,styles(mi),'Color',colors(mi,:),'DisplayName',labels(mi));
        end
        xlabel(xtext);ylabel('平均Oracle损失 (dB)');grid on;legend('Location','southoutside','NumColumns',2);
        title(array+"：L=32，包含先验违反");
    end
    finish(f,array+"_boundaries");
end
f=figure('Visible','off','Position',[100,100,1600,1200]);tiledlayout(2,2);
for array=g.arrayTypes
    for L=[32,1024]
        nexttile;hold on;
        for mi=1:numel(labels)
            q=T(T.arrayType==array&T.block=="snapshots"&T.nominalPriorValid&T.L==L&T.method==g.methodNames(mi),:);
            x=sort(q.medianLossDb);y=(1:numel(x))/numel(x);
            stairs(x,y,'Color',colors(mi,:),'LineStyle',extractBefore(styles(mi),strlength(styles(mi))), ...
                'DisplayName',labels(mi));
        end
        xlabel('几何中位Oracle损失 (dB)');ylabel('几何经验累积分布');
        title(sprintf('%s：L=%d，名义先验有效',array,L));grid on;legend('Location','southoutside','NumColumns',2);
    end
end
finish(f,"geometry_loss_cdf");
P=analysis.primary;
if width(P)>0
    f=figure('Visible','off','Position',[100,100,1300,850]);tiledlayout(2,2);
    for array=g.arrayTypes
        for other=["old_static","grid_covfit"]
            nexttile;q=P(P.arrayType==array&P.methodA=="raw_rec"&P.methodB==other,:);
            errorbar(q.L,q.differenceLossDb,q.differenceLossDb-q.simultaneousLowDb, ...
                q.simultaneousHighDb-q.differenceLossDb,'-o');yline(0,'k--');
            set(gca,'XScale','log');xlabel('快拍数L');ylabel('损失差 (dB)');
            xlim([min(q.L)/1.15,max(q.L)*1.15]);
            limits=[q.simultaneousLowDb;q.simultaneousHighDb;0];limits=limits(isfinite(limits));
            span=max(limits)-min(limits);if span>0,ylim([min(limits)-.1*span,max(limits)+.1*span]);end
            title([array+"：Rec 对 "+labels(g.methodNames==other);"族内同时95%区间"]);grid on;
        end
    end
    finish(f,"primary_simultaneous_intervals");
end
    function finish(f,name)
        set(findall(f,'-property','FontSize'),'FontSize',manifest.cfg.plot.fontSize);
        set(findall(f,'Type','line'),'LineWidth',manifest.cfg.plot.lineWidth);
        set(findall(f,'Type','errorbar'),'LineWidth',manifest.cfg.plot.lineWidth);
        savefig(f,fullfile(folder,name+".fig"));
        exportgraphics(f,fullfile(folder,name+".png"),'Resolution',manifest.cfg.plot.resolutionDpi);
        exportgraphics(f,fullfile(folder,name+".pdf"),'ContentType','vector');close(f);
    end
end
