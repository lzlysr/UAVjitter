function value = substream(cfg, experimentIndex, pointIndex, trialIndex, periodIndex, kind)
%SUBSTREAM Allocate one deterministic, nonoverlapping experiment substream.

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
value = experimentIndex * scheme.experimentStride + ...
    pointIndex * scheme.pointStride + trialIndex * scheme.trialStride + offset;
end
