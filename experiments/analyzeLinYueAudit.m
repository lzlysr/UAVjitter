function sensitivity = analyzeLinYueAudit(matPath)
%ANALYZELINYUEAUDIT 对NNLS审计异常做透明的中位数敏感性核查。
%   主结果保持原样；这里只对照“全部返回值”和“去掉KKT异常返回值”。
%   删除异常的结果不代替正确求解，也不能当作修复后的算法结果。
loaded=load(matPath,'result');r=loaded.result;cfg=r.cfg;
prefix="lin_yue_"+cfg.profile+"_nnls_audit";
loaded=load(fullfile(cfg.output.rawDir,prefix+".mat"),'audit');audit=loaded.audit;
detail=audit.detail;bad=detail(detail.relativeKktResidual>audit.tolerance,:);
sensitivity=table('Size',[0,7], ...
    'VariableTypes',{'double','double','string','double','double','double','double'}, ...
    'VariableNames',{'jitterStdDeg','snapshotCount','method','invalidKktCount', ...
    'originalMedianDb','withoutInvalidMedianDb','sensitivityDb'});
rows=cell(height(bad),1);
for q=1:height(bad)
    ji=find(r.jitterStdDeg==bad.jitterStdDeg(q));
    ip=find(r.scanPoints==bad.snapshotCount(q));
    mi=find(r.methodNames==bad.method(q));
    values=r.raw.sinrLinear(:,:,ip,ji,mi);mask=true(size(values));
    same=bad.jitterStdDeg==bad.jitterStdDeg(q) & bad.snapshotCount==bad.snapshotCount(q) & bad.method==bad.method(q);
    mask(sub2ind(size(values),bad.trial(same),bad.period(same)))=false;
    allMedian=10*log10(median(values(:)));
    validMedian=10*log10(median(values(mask)));
    rows{q}=table(bad.jitterStdDeg(q),bad.snapshotCount(q),bad.method(q), ...
        sum(same),allMedian,validMedian,validMedian-allMedian, ...
        'VariableNames',{'jitterStdDeg','snapshotCount','method','invalidKktCount', ...
        'originalMedianDb','withoutInvalidMedianDb','sensitivityDb'});
end
if ~isempty(rows),sensitivity=unique(vertcat(rows{:}));end
writetable(sensitivity,fullfile(cfg.output.tableDir,prefix+"_sensitivity.csv"));
disp(sensitivity);
end
