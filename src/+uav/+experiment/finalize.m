function result = finalize(result, summaryTable, fig, cfg, actualCounts, timerValue)
%FINALIZE Attach metadata, export artifacts, close the figure, and return paths.

runtimeSec = toc(timerValue);
result.cfg = cfg;
result.metadata = uav.experiment.metadata(cfg, result.experimentName, ...
    actualCounts, runtimeSec);
result.timestamp = result.metadata.generatedAt;
result.runtimeSec = runtimeSec;
result.status = "complete";
artifacts = uav.experiment.writeArtifacts(result, summaryTable, fig, cfg);
artifacts.status = "complete";
result.artifacts = artifacts;
close(fig);
end
