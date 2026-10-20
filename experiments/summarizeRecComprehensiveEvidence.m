function aggregates=summarizeRecComprehensiveEvidence(directory)
%SUMMARIZERECCOMPREHENSIVEEVIDENCE 等权汇总几何，不把试次数膨胀成几何数。
%   所有分层都导出；本表是描述性摘要，同时区间应读取primary表。
load(fullfile(directory,'analysis.mat'),'analysis');T=analysis.summary;
fields={'arrayType','block','geometryFamily','method','L','snrDb','inrDb','secondInrDb', ...
    'J','halfWidthDeg','jitterStdDeg','gridSize','population','nominalPriorValid'};
keys=strings(height(T),1);
for index=1:numel(fields),keys=keys+"|"+string(T.(fields{index}));end
[uniqueKeys,~,group]=unique(keys); rows=cell(numel(uniqueKeys),1);
for index=1:numel(uniqueKeys)
    q=T(group==index,:); row=q(1,fields);
    row.geometryCount=height(q);
    row.meanGeometryMedianLossDb=mean(q.medianLossDb);
    row.medianGeometryLossDb=median(q.medianLossDb);
    row.p90GeometryLossDb=quantile(q.medianLossDb,.9);
    row.maxGeometryLossDb=max(q.medianLossDb);
    row.meanGeometryMedianSinrDb=mean(q.medianSinrDb);
    row.meanOracleMedianDb=mean(q.oracleMedianDb);
    row.meanTrueIncmNominalMedianDb=mean(q.trueIncmNominalMedianDb);
    row.meanTrialOutageRate=mean(q.outageRate);
    finiteLoss=q.medianLossDb(isfinite(q.medianLossDb));
    if isempty(finiteLoss),row.geometryMedianLossOver3Rate=NaN;
    else,row.geometryMedianLossOver3Rate=mean(finiteLoss>3);end
    row.meanIncmNmse=mean(q.medianIncmNmse);
    row.meanShapeNmse=mean(q.medianShapeNmse);
    row.numericalFailures=sum(q.numericalFailures);row.successes=sum(q.successes);
    row.notApplicable=sum(q.notApplicable);row.kktViolations=sum(q.kktViolations);
    row.meanDesiredOutsideRate=mean(q.desiredOutsideRate);
    row.meanInterfererInsideRate=mean(q.interfererInsideRate);
    rows{index}=row;
end
aggregates=vertcat(rows{:});writetable(aggregates,fullfile(directory,'geometry_aggregates.csv'));
aggregationSource=fileread(mfilename('fullpath')+".m");
save(fullfile(directory,'geometry_aggregates.mat'),'aggregates','aggregationSource');
% 同世界几何的阵型配对表：两阵型的有效先验子集不同，不能直接比较不同总体的均值。
base=T(T.block=="snapshots",{'geometryId','L','method','nominalPriorValid', ...
    'medianLossDb','oracleMedianDb','medianSinrDb','arrayType'});
if ~isempty(base)
    first=removevars(base(base.arrayType=="ULA",:),'arrayType');
    second=removevars(base(base.arrayType=="UPA",:),'arrayType');
    first.Properties.VariableNames(4:end)=cellstr(string(first.Properties.VariableNames(4:end))+"_ULA");
    second.Properties.VariableNames(4:end)=cellstr(string(second.Properties.VariableNames(4:end))+"_UPA");
    matched=innerjoin(first,second,'Keys',{'geometryId','L','method'});
    matched.bothPriorsValid=matched.nominalPriorValid_ULA & matched.nominalPriorValid_UPA;
    matched.lossUPAMinusULADb=matched.medianLossDb_UPA-matched.medianLossDb_ULA;
    writetable(matched,fullfile(directory,'array_matched_geometry.csv'));
end
end
