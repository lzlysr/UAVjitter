function publishStagedBundle(staged, artifacts, backups, operations)
%PUBLISHSTAGEDBUNDLE Replace a result bundle with a commit point and rollback.

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
if nargin < 4
    operations = struct('move', @moveChecked, 'delete', @deleteChecked);
elseif ~isstruct(operations) || ~all(isfield(operations, {'move', 'delete'})) || ...
        ~isa(operations.move, 'function_handle') || ...
        ~isa(operations.delete, 'function_handle')
    error('uav:experiment:InvalidArtifactOperations', ...
        'Publication operations must provide move and delete function handles.');
end

backedUp = false(size(fields));
published = false(size(fields));
try
    for index = 1:numel(fields)
        field = fields{index};
        if isfile(artifacts.(field))
            operations.move(artifacts.(field), backups.(field));
            backedUp(index) = true;
        end
    end
    for index = 1:numel(fields)
        field = fields{index};
        operations.move(staged.(field), artifacts.(field));
        published(index) = true;
    end
catch publishException
    rollbackFailures = rollbackBeforeCommit(fields, staged, artifacts, ...
        backups, backedUp, published, operations);
    if isempty(rollbackFailures)
        rethrow(publishException);
    end
    message = sprintf(['Artifact publication failed and rollback was incomplete. ' ...
        'Preserved recovery files: %s'], strjoin(rollbackFailures, ', '));
    rollbackException = MException( ...
        'uav:experiment:ArtifactRollbackFailure', '%s', message);
    rollbackException = addCause(rollbackException, publishException);
    throw(rollbackException);
end

% Commit point: every new final exists. Cleanup can no longer trigger rollback.
for index = 1:numel(fields)
    field = fields{index};
    if backedUp(index) && isfile(backups.(field))
        try
            operations.delete(backups.(field));
        catch cleanupException
            warning('uav:experiment:BackupCleanupFailure', ...
                'Published bundle is complete; retained backup %s (%s).', ...
                backups.(field), cleanupException.message);
        end
    end
end
end

function failures = rollbackBeforeCommit(fields, staged, artifacts, backups, ...
        backedUp, published, operations)
failures = strings(0, 1);
finalCleared = true(size(fields));
for index = 1:numel(fields)
    field = fields{index};
    if published(index) && isfile(artifacts.(field))
        try
            operations.delete(artifacts.(field));
        catch
            finalCleared(index) = false;
            failures(end + 1, 1) = string(backups.(field)); %#ok<AGROW>
        end
    end
end
for index = 1:numel(fields)
    field = fields{index};
    if backedUp(index) && isfile(backups.(field))
        if finalCleared(index) && ~isfile(artifacts.(field))
            try
                operations.move(backups.(field), artifacts.(field));
            catch
                failures(end + 1, 1) = string(backups.(field)); %#ok<AGROW>
            end
        elseif ~any(failures == string(backups.(field)))
            failures(end + 1, 1) = string(backups.(field)); %#ok<AGROW>
        end
    end
end
for index = 1:numel(fields)
    field = fields{index};
    if isfile(staged.(field))
        try
            operations.delete(staged.(field));
        catch
            failures(end + 1, 1) = string(staged.(field)); %#ok<AGROW>
        end
    end
end
failures = unique(failures, 'stable');
end

function moveChecked(source, destination)
[moved, message] = movefile(source, destination, 'f');
if ~moved
    error('uav:experiment:ArtifactPublishFailure', ...
        'Could not publish %s: %s', destination, message);
end
end

function deleteChecked(pathValue)
try
    delete(pathValue);
catch exception
    cleanupException = MException( ...
        'uav:experiment:ArtifactCleanupFailure', ...
        'Could not delete %s.', pathValue);
    cleanupException = addCause(cleanupException, exception);
    throw(cleanupException);
end
end
