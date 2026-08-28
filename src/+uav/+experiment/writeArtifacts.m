function artifacts = writeArtifacts(result, summaryTable, fig, cfg)
%WRITEARTIFACTS 暂存并发布一套完整、可复现的实验结果文件。
%   六种格式先写入各目标目录中的临时文件；任何导出失败都会清理暂存
%   文件，不会覆盖上一轮完整结果。全部导出成功后再交给事务发布器。

if ~isstruct(result) || ~isfield(result, 'metadata') || ...
        ~istable(summaryTable) || ~isgraphics(fig, 'figure')
    error('uav:experiment:InvalidArtifactInput', ...
        'Result, summary table, and figure are required.');
end
uav.util.ensureDirectories(cfg.output);
baseName = result.metadata.experimentName + "_" + cfg.profile;
artifacts = finalPaths(baseName, cfg);
finalValues = string(struct2cell(artifacts));
if ~cfg.output.overwrite && any(isfile(finalValues))
    error('uav:experiment:ArtifactExists', ...
        'Configured outputs already exist and overwrite is disabled.');
end

staged = stagingPaths(artifacts);
backups = stagingPaths(artifacts);
% 外层只负责暂存文件；备份的生命周期由publishStagedBundle管理。
stageCleanup = onCleanup(@() deleteExisting(string(struct2cell(staged))));
result.status = "complete";
result.artifacts = artifacts;
% MAT保留完整原始数组；CSV用于表格；FIG/PDF/PNG满足编辑、印刷和预览。
save(staged.mat, 'result', '-v7.3');
writetable(summaryTable, staged.csv);
savefig(fig, staged.fig);
exportgraphics(fig, staged.pdf, 'ContentType', 'vector');
exportgraphics(fig, staged.png, 'Resolution', cfg.plot.resolutionDpi);
manifest = struct( ...
    'status', "complete", ...
    'experimentName', result.metadata.experimentName, ...
    'profile', cfg.profile, ...
    'generatedAt', result.metadata.generatedAt, ...
    'artifacts', artifacts);
writeJson(staged.log, manifest);
% 只有六个暂存文件都存在时才进入最终发布阶段。
uav.experiment.publishStagedBundle(staged, artifacts, backups);
end

function artifacts = finalPaths(baseName, cfg)
artifacts = struct( ...
    'mat', string(fullfile(cfg.output.rawDir, baseName + ".mat")), ...
    'csv', string(fullfile(cfg.output.tableDir, baseName + ".csv")), ...
    'fig', string(fullfile(cfg.output.figureDir, baseName + ".fig")), ...
    'pdf', string(fullfile(cfg.output.figureDir, baseName + ".pdf")), ...
    'png', string(fullfile(cfg.output.figureDir, baseName + ".png")), ...
    'log', string(fullfile(cfg.output.logDir, baseName + ".json")));
end

function staged = stagingPaths(artifacts)
% 临时文件与最终文件位于同一目录，使最终movefile保持同卷、快速且可控。
fields = fieldnames(artifacts);
staged = struct();
for index = 1:numel(fields)
    field = fields{index};
    finalPath = artifacts.(field);
    [directory, ~, extension] = fileparts(finalPath);
    staged.(field) = string(tempname(directory)) + extension;
end
end

function writeJson(pathValue, manifest)
fileId = fopen(pathValue, 'w');
if fileId < 0
    error('uav:experiment:LogWriteFailure', 'Could not create the manifest log.');
end
cleanup = onCleanup(@() fclose(fileId));
fwrite(fileId, jsonencode(manifest, 'PrettyPrint', true), 'char');
end

function deleteExisting(paths)
for index = 1:numel(paths)
    if isfile(paths(index))
        delete(paths(index));
    end
end
end
