function artifacts = writeArtifacts(result, summaryTable, fig, cfg)
%WRITEARTIFACTS Stage and publish one complete deterministic result bundle.

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
stageCleanup = onCleanup(@() deleteExisting([string(struct2cell(staged)); ...
    string(struct2cell(backups))]));
result.status = "complete";
result.artifacts = artifacts;
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
publishBundle(staged, artifacts, backups);
end

function publishBundle(staged, artifacts, backups)
fields = fieldnames(artifacts);
backedUp = false(size(fields));
published = false(size(fields));
try
    for index = 1:numel(fields)
        field = fields{index};
        if isfile(artifacts.(field))
            moveChecked(artifacts.(field), backups.(field));
            backedUp(index) = true;
        end
    end
    for index = 1:numel(fields)
        field = fields{index};
        moveChecked(staged.(field), artifacts.(field));
        published(index) = true;
    end
    deleteExisting(string(struct2cell(backups)));
catch exception
    for index = 1:numel(fields)
        field = fields{index};
        if published(index) && isfile(artifacts.(field))
            delete(artifacts.(field));
        end
    end
    for index = 1:numel(fields)
        field = fields{index};
        if backedUp(index) && isfile(backups.(field))
            movefile(backups.(field), artifacts.(field), 'f');
        end
    end
    rethrow(exception);
end
end

function moveChecked(source, destination)
[moved, message] = movefile(source, destination, 'f');
if ~moved
    error('uav:experiment:ArtifactPublishFailure', ...
        'Could not publish %s: %s', destination, message);
end
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
