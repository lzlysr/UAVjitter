function publishStagedBundle(staged, artifacts, backups)
%PUBLISHSTAGEDBUNDLE Atomically replace a result bundle with rollback.

if ~isstruct(staged) || ~isstruct(artifacts) || ~isstruct(backups) || ...
        ~isequal(fieldnames(staged), fieldnames(artifacts), fieldnames(backups))
    error('uav:experiment:InvalidArtifactBundle', ...
        'Staged, final, and backup bundles must have identical fields.');
end
fields = fieldnames(artifacts);
if isempty(fields) || any(~isfile(string(struct2cell(staged))))
    error('uav:experiment:InvalidArtifactBundle', ...
        'Every staged artifact must exist before publication.');
end

backedUp = false(size(fields));
published = false(size(fields));
cleanup = onCleanup(@() deleteExisting([string(struct2cell(staged)); ...
    string(struct2cell(backups))]));
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
            moveChecked(backups.(field), artifacts.(field));
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

function deleteExisting(paths)
for index = 1:numel(paths)
    if isfile(paths(index))
        delete(paths(index));
    end
end
end
