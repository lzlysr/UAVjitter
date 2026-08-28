function value = substream(cfg, experimentIndex, pointIndex, trialIndex, periodIndex, kind)
%SUBSTREAM 为实验层级分配确定性且互不重叠的随机子流编号。
%   编号层次为：实验 -> 扫描点 -> 试次 -> 用途/周期。stride在config中
%   经过范围校验，确保scenario、snapshot和bootstrap子流不会碰撞。

indices = [experimentIndex, pointIndex, trialIndex, periodIndex];
if any(~isfinite(indices)) || any(indices < 0) || any(indices ~= floor(indices))
    error('uav:experiment:InvalidSubstreamIndex', ...
        'Experiment substream indices must be finite nonnegative integers.');
end
scheme = cfg.experiments.substreams;
kind = lower(string(kind));
switch kind
    case "scenario"
        offset = scheme.scenarioOffset;
    case "snapshot"
        offset = scheme.snapshotOffset + periodIndex;
    case "bootstrap"
        offset = scheme.bootstrapOffset;
    otherwise
        error('uav:experiment:InvalidSubstreamKind', ...
            'Substream kind must be scenario, snapshot, or bootstrap.');
end
% 这里只计算编号，不创建RandStream，也不改变全局rng。
value = experimentIndex * scheme.experimentStride + ...
    pointIndex * scheme.pointStride + trialIndex * scheme.trialStride + offset;
end
