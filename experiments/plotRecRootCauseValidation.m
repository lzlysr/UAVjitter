function paths = plotRecRootCauseValidation(resultFile, keepOpen, selectedNames)
%PLOTRECROOTCAUSEVALIDATION 从既有验证原始数据绘制诊断图，不重新仿真。
%   第一组图显示试次周期线性均值的中位数与试次级bootstrap 95%区间；
%   第二组图显示全部周期的5%分位；第三组图显示低于15 dB的观测比例。
%   后两组只画样本统计量，不将它们冒充置信区间。
%   默认只展示鲁棒基线与两种细化候选；全部11种方法始终保存在MAT/CSV，
%   可通过selectedNames显式选择其他方法。keepOpen=false用于批处理关闭图窗。
arguments
    resultFile (1,1) string
    keepOpen (1,1) logical = true
    selectedNames (1,:) string = ["oracle", "esp", "wco", "dl_mvdr", ...
        "upa_covfit_incm", "refined_candidate_covfit", "refined_narrow_sector"]
end
saved = load(resultFile, 'results');
r = saved.results;
cfg = r.cfg;
study = cfg.experiments.recRootCause;
[present, indices] = ismember(selectedNames, r.methodNames);
assert(all(present) && numel(unique(selectedNames)) == numel(selectedNames), ...
    '绘图方法必须属于结果中已有的唯一方法集合。');
directory = fullfile(fileparts(resultFile), "figures_" + cfg.profile);
if ~isfolder(directory), mkdir(directory); end
paths = strings(0, 1);
stream = uav.util.makeStream(cfg.random.baseSeed, study.validationSeedOffset*3);
cfg.plot.visible = keepOpen;
catalogNames = ["oracle", "esp", "wco", "gu_leshem2012_rec", "upa_covfit_incm", ...
    "candidate_covfit", "dl_mvdr", "old_static", "true_incm_nominal", ...
    "refined_candidate_covfit", "refined_narrow_sector"];
catalogLabels = ["理想上界", "ESP", "WCO", "原式Rec", "原CovFit", ...
    "未细化候选", "DL", "旧static", "真实INCM＋名义方向", "细化CovFit（25°）", "细化CovFit（7°）"];
styles = ["--", ":", "-", "--", "-.", ":", "-.", "--", ":", "-", "-"];
markers = ["none", "p", "h", "x", "*", "d", "+", "^", "v", "s", "o"];
colors = lines(numel(catalogNames));
colors(1, :) = zeros(1, 3);
for snapshotIndex = 1:numel(study.validationSnapshotCounts)
    L = study.validationSnapshotCounts(snapshotIndex);
    center = nan(numel(study.validationJitterDeg), numel(indices));
    lower = center; upper = center; low = center; outage = center;
    for point = 1:numel(study.validationJitterDeg)
        for method = 1:numel(indices)
            raw = r.sinrLinear(:, :, point, snapshotIndex, indices(method));
            if all(isnan(raw), 'all'), continue; end
            trialMean = mean(raw, 2);
            resample = randi(stream, size(raw, 1), size(raw, 1), cfg.experiments.bootstrapReplicates);
            bootMedian = median(trialMean(resample), 1);
            interval = prctile(bootMedian, 100*cfg.experiments.aggregation.confidenceQuantiles);
            center(point, method) = 10*log10(median(trialMean));
            lower(point, method) = 10*log10(interval(1));
            upper(point, method) = 10*log10(interval(2));
            low(point, method) = 10*log10(prctile(raw(:), 100*cfg.experiments.aggregation.lowQuantile));
            outage(point, method) = mean(raw < 10^(study.validationOutageThresholdDb/10), 'all');
        end
    end
    for panel = ["median", "period_low", "outage"]
        fig = uav.plot.newFigure(cfg);
        ax = axes(fig); hold(ax, 'on');
        labels = strings(size(selectedNames));
        for method = 1:numel(indices)
            location = find(catalogNames == selectedNames(method));
            labels(method) = catalogLabels(location);
            if selectedNames(method) == "gu_leshem2012_rec" && L < cfg.derived.M
                labels(method) = labels(method) + "（L<M不适用）";
            end
            if panel == "median"
                errorbar(ax, study.validationJitterDeg, center(:, method), ...
                    center(:, method)-lower(:, method), upper(:, method)-center(:, method), ...
                    'LineStyle', styles(location), 'Marker', markers(location), ...
                    'Color', colors(location, :), 'LineWidth', cfg.plot.lineWidth, ...
                    'MarkerSize', cfg.plot.markerSize);
            else
                if panel == "period_low", values = low(:, method); else, values = outage(:, method); end
                plot(ax, study.validationJitterDeg, values, 'LineStyle', styles(location), ...
                    'Marker', markers(location), 'Color', colors(location, :), ...
                    'LineWidth', cfg.plot.lineWidth, 'MarkerSize', cfg.plot.markerSize);
            end
        end
        xlabel(ax, '三轴姿态抖动RMS（°）', 'FontSize', cfg.plot.fontSize);
        if panel == "median"
            ylabel(ax, '输出SINR中位数（dB）', 'FontSize', cfg.plot.fontSize);
        elseif panel == "period_low"
            ylabel(ax, '逐周期SINR第5百分位（dB）', 'FontSize', cfg.plot.fontSize);
        else
            ylabel(ax, sprintf('SINR低于%g dB的观测比例', study.validationOutageThresholdDb), ...
                'FontSize', cfg.plot.fontSize);
        end
        title(ax, sprintf('%s：L=%d，%d试次×%d周期', cfg.profile, L, ...
            study.validationTrials, study.validationPeriods), ...
            'FontSize', cfg.plot.fontSize, 'Interpreter', 'none');
        set(ax, 'FontSize', cfg.plot.fontSize); grid(ax, 'on'); box(ax, 'on');
        legend(ax, labels, 'Location', 'bestoutside', 'FontSize', cfg.plot.fontSize);
        stem = fullfile(directory, "validation_L" + L + "_" + panel);
        savefig(fig, stem + ".fig");
        exportgraphics(fig, stem + ".pdf", 'ContentType', 'vector');
        exportgraphics(fig, stem + ".png", 'Resolution', cfg.plot.resolutionDpi);
        paths = [paths; stem + [".fig"; ".pdf"; ".png"]]; %#ok<AGROW>
        if ~keepOpen, close(fig); end
    end
    % 取该L下目标方法最差周期所处的RMS，展示该场景的全部周期CDF。
    % 用对数概率轴显示1%以下尾部，避免中位数和5%分位掩盖极端失败。
    targetIndex = find(r.methodNames=="refined_narrow_sector");
    targetValues = r.sinrLinear(:,:,:,snapshotIndex,targetIndex);
    [~,worstIndex] = min(targetValues(:));
    [~,~,worstPoint] = ind2sub(size(targetValues),worstIndex);
    fig = uav.plot.newFigure(cfg);
    ax = axes(fig); hold(ax,'on');
    labels = strings(size(selectedNames));
    for method = 1:numel(indices)
        location = find(catalogNames==selectedNames(method));
        labels(method) = catalogLabels(location);
        values = r.sinrLinear(:,:,worstPoint,snapshotIndex,indices(method));
        if all(isnan(values),'all')
            labels(method) = labels(method)+"（不适用）";
        end
        [x,p] = uav.plot.finiteCdf(10*log10(values(:)));
        % Stair对象不支持MarkerIndices，先取阶梯坐标再用Line绘制同一经验CDF。
        if isempty(x), stepX = x; stepP = p; else, [stepX,stepP] = stairs(x,p); end
        markerIndices = unique(round(linspace(1,max(numel(stepX),1),numel(study.validationJitterDeg))));
        if isempty(x), markerIndices = []; end
        plot(ax,stepX,stepP,'LineStyle',styles(location),'Marker',markers(location), ...
            'MarkerIndices',markerIndices,'Color',colors(location,:), ...
            'LineWidth',cfg.plot.lineWidth,'MarkerSize',cfg.plot.markerSize);
    end
    set(ax,'YScale','log','FontSize',cfg.plot.fontSize); grid(ax,'on'); box(ax,'on');
    xlabel(ax,'逐周期输出SINR（dB）','FontSize',cfg.plot.fontSize);
    ylabel(ax,'经验累积概率（对数纵轴）','FontSize',cfg.plot.fontSize);
    title(ax,sprintf('L=%d，最差样本对应RMS=%g°',L,study.validationJitterDeg(worstPoint)), ...
        'FontSize',cfg.plot.fontSize);
    legend(ax,labels,'Location','bestoutside','FontSize',cfg.plot.fontSize);
    stem = fullfile(directory,"validation_L"+L+"_worst_case_cdf");
    savefig(fig,stem+".fig");
    exportgraphics(fig,stem+".pdf",'ContentType','vector');
    exportgraphics(fig,stem+".png",'Resolution',cfg.plot.resolutionDpi);
    paths = [paths;stem+[".fig";".pdf";".png"]]; %#ok<AGROW>
    if ~keepOpen, close(fig); end
end
end
