function [RinHat, nextState, diagnostics] = temporalCovfit(obs, priorState, cfg, baseIncm, base, method)
%TEMPORALCOVFIT 跟踪干扰候选并在当前字典上施加非负功率连续性约束。
%   首周期沿用细化CovFit。后续用过去已拟合的干扰确认当前SCM提取的候选。
%   最终期望导向固定，历史按试次独立初始化。

if nargin < 6, method = "incm_covfit_refined_temporal"; end
if nargin < 4 || isempty(baseIncm)
    [baseIncm, base] = uav.incm.localCovfit(obs,cfg);
end
RinHat = baseIncm;
diagnostics = base;
diagnostics.method = method;
diagnostics.temporalApplied = false;
nodes = base.interfererCandidateNodes;
powers = base.interfererFittedPowers;
nextState = makeState(nodes, powers, cfg);
if isempty(priorState) || isempty(priorState.nodes) || min(cfg.array.Nx,cfg.array.Ny) < 2
    return;
end

settings = cfg.incm.temporalCovfit;
switch method
    case "incm_covfit_refined_temporal_capon"
        [proposed, proposedPowers] = uav.spectrum.fullCaponNodes(obs,base.noisePower,cfg);
    case "incm_covfit_refined_temporal_hybrid"
        [proposed, proposedPowers] = uav.spectrum.hybridCaponNodes(obs,base.noisePower,cfg);
    otherwise
        [proposed, proposedPowers] = uav.spectrum.shiftInvariantNodes(obs.Rx,base.noisePower,cfg);
end
diagnostics.proposedCandidateNodes = proposed;
pairs = associate(priorState, proposed, proposedPowers, settings);
if isempty(pairs)
    return;
end
desired = base.desiredCandidateNodes;
current = proposed(pairs(:,2),:);
transported = priorState.powers(pairs(:,1));
% 走到固定排除区内的候选须由历史确认，且仍要与当前期望字典分开；
% 名义期望导向从未参与追踪或更新。
duplicateDesired = any(distances(current, desired) <= settings.mergeDistance,2);
current(duplicateDesired,:) = [];
transported(duplicateDesired) = [];
pairs(duplicateDesired,:) = [];
if isempty(current)
    return;
end

% 保留相对强且未被关联结果覆盖的当前CovFit候选，容许新干扰进入。
newStrong = powers > settings.minimumTrackPowerRatio * max([powers;0]);
notCovered = all(distances(nodes,current) > settings.associationDistance,2);
fresh = nodes(newStrong & notCovered,:);
room = max(0,cfg.incm.localCovfit.maximumCandidateCount-size(current,1));
fresh = fresh(1:min(size(fresh,1),room),:);
nodes = [current;fresh];
transported = [transported;zeros(size(fresh,1),1)];
matched = [true(size(current,1),1);false(size(fresh,1),1)];

% 只对成功关联的干扰系数加入二次连续性项。期望候选联合拟合后剔除。
D = size(desired,1);
A = uav.array.steeringUPA([desired;nodes].',cfg.array.Nx,cfg.array.Ny, ...
    cfg.array.spacingOverLambda);
penalty = cfg.derived.M^2 * settings.powerContinuityWeight * [false(D,1);matched];
target = [zeros(D,1);transported];
[p, fit] = uav.incm.fitNonnegativePowers(obs.Rx,A,base.noisePower, ...
    cfg.incm.localCovfit.fitRegularization,penalty,target);
kkt = uav.incm.nonnegativeFitKkt(obs.Rx,A,base.noisePower,p, ...
    cfg.incm.localCovfit.fitRegularization,penalty,target);
if kkt > cfg.incm.localCovfit.kktTolerance
    error('uav:incm:InvalidTemporalCovfitSolution','时序CovFit功率解未满足非负KKT条件。');
end
Ai = A(:,D+1:end);
powers = p(D+1:end);
RinHat = (Ai.*powers.')*Ai' + base.noisePower*eye(cfg.derived.M);
RinHat = (RinHat + RinHat')/2;
nextState = makeState(nodes, powers, cfg);
diagnostics.candidateManifold = A;
diagnostics.interfererCandidateNodes = nodes;
diagnostics.fittedPowers = p;
diagnostics.desiredFittedPowers = p(1:D);
diagnostics.interfererFittedPowers = powers;
diagnostics.fitDiagnostics = fit;
diagnostics.relativeKktResidual = kkt;
diagnostics.temporalApplied = true;
diagnostics.temporalPairs = pairs;
diagnostics.temporalMatched = matched;
diagnostics.transportedPowers = transported;
end

function state = makeState(nodes, powers, cfg)
% 只传递本周期拟合可靠的主干扰，避免弱伪峰在状态中长期存在。
reference = max([powers;0]);
keep = powers > cfg.incm.temporalCovfit.minimumTrackPowerRatio * reference;
state = struct('nodes', nodes(keep,:), 'powers', powers(keep));
end

function pairs = associate(history, current, currentPowers, settings)
% 方向余弦和功率都须与历史相容，再按最邻近距离执行一对一分配。
cost = distances(history.nodes,current);
ratio = currentPowers.' ./ history.powers;
cost(ratio < 1/settings.maximumPowerRatio | ratio > settings.maximumPowerRatio) = inf;
pairs = zeros(0,2);
while ~isempty(cost)
    [closest, index] = min(cost, [], 'all', 'linear');
    if closest > settings.associationDistance
        break;
    end
    [previous, present] = ind2sub(size(cost), index);
    pairs(end+1,:) = [previous, present]; %#ok<AGROW>
    cost(previous,:) = inf;
    cost(:,present) = inf;
end
end

function distance = distances(left,right)
distance = hypot(left(:,1)-right(:,1).',left(:,2)-right(:,2).');
end
