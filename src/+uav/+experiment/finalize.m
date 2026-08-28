function result = finalize(result, summaryTable, fig, cfg, actualCounts, timerValue)
%FINALIZE 补充元数据、导出完整结果包并关闭图窗。
%   onCleanup在函数刚开始就注册，因此writeArtifacts报错时也会关闭图窗，
%   防止批量实验因隐藏图窗累计而耗尽系统资源。

runtimeSec = toc(timerValue);
figureCleanup = onCleanup(@() closeIfValid(fig));
% 保存的是实验实际使用的完整cfg，其中包含固定标定后的epsilonS。
result.cfg = cfg;
result.metadata = uav.experiment.metadata(cfg, result.experimentName, ...
    actualCounts, runtimeSec);
result.timestamp = result.metadata.generatedAt;
result.runtimeSec = runtimeSec;
artifacts = uav.experiment.writeArtifacts(result, summaryTable, fig, cfg);
% 只有全部文件成功发布后，返回结果才标记为complete。
artifacts.status = "complete";
result.status = "complete";
result.artifacts = artifacts;
end

function closeIfValid(fig)
if isgraphics(fig, 'figure')
    close(fig);
end
end
