function result = finalize(result, summaryTable, fig, cfg, actualCounts, timerValue)
%FINALIZE Attach metadata, export artifacts, close the figure, and return paths.

runtimeSec = toc(timerValue);
figureCleanup = onCleanup(@() closeIfValid(fig));
result.cfg = cfg;
result.metadata = uav.experiment.metadata(cfg, result.experimentName, ...
    actualCounts, runtimeSec);
result.timestamp = result.metadata.generatedAt;
result.runtimeSec = runtimeSec;
artifacts = uav.experiment.writeArtifacts(result, summaryTable, fig, cfg);
artifacts.status = "complete";
result.status = "complete";
result.artifacts = artifacts;
end

function closeIfValid(fig)
if isgraphics(fig, 'figure')
    close(fig);
end
end
