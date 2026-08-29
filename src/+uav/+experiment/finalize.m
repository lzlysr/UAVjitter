function result = finalize(result, summaryTable, fig, cfg, actualCounts, timerValue)
%FINALIZE 补充元数据、导出结果包，并按运行方式管理图窗生命周期。
%   导出失败时始终关闭图窗；导出成功后由cfg.plot.closeAfterExport决定
%   是否关闭。单实验默认保留，run_all批处理会把该开关设为true。

runtimeSec = toc(timerValue);
try
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
catch exception
    % 失败图窗不具有完整实验结果，始终关闭后再原样抛出异常。
    closeIfValid(fig);
    rethrow(exception);
end

if cfg.plot.closeAfterExport
    closeIfValid(fig);
end
end

function closeIfValid(fig)
if isgraphics(fig, 'figure')
    close(fig);
end
end
