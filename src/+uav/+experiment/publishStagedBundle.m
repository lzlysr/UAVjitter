function publishStagedBundle(staged, artifacts, backups, operations)
%PUBLISHSTAGEDBUNDLE 带明确提交点地替换一整套实验结果。
%   staged/artifacts/backups具有完全相同字段，分别表示暂存文件、最终文件
%   和旧结果备份。正常调用使用真实move/delete操作；可选operations只为
%   故障注入测试提供替代操作。
%
%   事务语义：
%   1) 先把所有旧final移动到backup；
%   2) 再把所有staged移动到final；
%   3) 第2步全部成功即到达提交点，此后绝不回滚完整新结果；
%   4) 提交点之前失败则删除已发布新文件并恢复旧备份。

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
    % 阶段一：保护所有已经存在的旧结果。
    for index = 1:numel(fields)
        field = fields{index};
        if isfile(artifacts.(field))
            operations.move(artifacts.(field), backups.(field));
            backedUp(index) = true;
        end
    end
    % 阶段二：逐个发布已完整生成的暂存文件。
    for index = 1:numel(fields)
        field = fields{index};
        operations.move(staged.(field), artifacts.(field));
        published(index) = true;
    end
catch publishException
    % 只有尚未到达提交点时才允许回滚。
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

% 提交点：所有新final均已存在。此后的备份清理失败不能触发回滚。
for index = 1:numel(fields)
    field = fields{index};
    if backedUp(index) && isfile(backups.(field))
        try
            operations.delete(backups.(field));
        catch cleanupException
            % 保留未删除备份并报警，完整的新结果仍然有效。
            warning('uav:experiment:BackupCleanupFailure', ...
                'Published bundle is complete; retained backup %s (%s).', ...
                backups.(field), cleanupException.message);
        end
    end
end
end

function failures = rollbackBeforeCommit(fields, staged, artifacts, backups, ...
        backedUp, published, operations)
%ROLLBACKBEFORECOMMIT 尽最大努力回到发布前状态，并保留无法恢复的备份。
failures = strings(0, 1);
finalCleared = true(size(fields));
% 先删除已经发布的新文件，为恢复同名旧文件腾出目标路径。
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
% 仅当目标路径已清空时恢复旧文件；恢复失败的backup绝不能删除。
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
% 清理尚未发布的暂存文件；清理失败路径同样报告给调用者。
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
