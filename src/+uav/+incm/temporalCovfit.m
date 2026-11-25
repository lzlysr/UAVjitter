function [RinHat, nextState, diagnostics] = temporalCovfit(obs, priorState, cfg, baseIncm, base, method)
%TEMPORALCOVFIT 跟踪干扰候选并在当前字典上施加非负功率连续性约束。
%   首周期用各自的单周期CovFit。后续用上期干扰确认当期方向候选。
%   最终期望导向固定，历史按试次独立初始化。

if nargin < 6, method = "incm_covfit_refined_temporal"; end
settings = cfg.incm.temporalCovfit;
pureSubspace = method == "incm_covfit_refined_temporal_pure_subspace";
resolveConflict = ismember(method,["incm_covfit_refined_temporal", ...
    "incm_covfit_refined_temporal_capon_base","incm_covfit_refined_temporal_pure_subspace"]) && ...
    isfield(settings,'resolveDesiredPeakConflict') && settings.resolveDesiredPeakConflict;
if pureSubspace
    % 纯子空间版独立生成底座，首周期和关联失败时也使用这一解。
    [baseIncm, base] = uav.incm.subspaceCovfit(obs,cfg);
elseif nargin < 4 || isempty(baseIncm)
    [baseIncm, base] = uav.incm.localCovfit(obs,cfg);
end
RinHat = baseIncm;
diagnostics = base;
diagnostics.method = method;
diagnostics.temporalApplied = false;
diagnostics.reclassifiedDesiredPeak = zeros(0,2);
nodes = base.interfererCandidateNodes;
powers = base.interfererFittedPowers;
nextState = makeState(nodes, powers, cfg, resolveConflict);
if isempty(priorState) || isempty(priorState.nodes) || min(cfg.array.Nx,cfg.array.Ny) < 2
    return;
end

switch method
    case {"incm_covfit_refined_temporal","incm_covfit_refined_temporal_capon"}
        [proposed, proposedPowers] = uav.spectrum.fullCaponNodes(obs,base.noisePower,cfg);
    case "incm_covfit_refined_temporal_hybrid"
        [proposed, proposedPowers] = uav.spectrum.hybridCaponNodes(obs,base.noisePower,cfg);
    case "incm_covfit_refined_temporal_capon_base"
        [proposed, proposedPowers] = uav.spectrum.shiftInvariantNodes(obs.Rx,base.noisePower,cfg);
    case "incm_covfit_refined_temporal_pure_subspace"
        proposed = base.proposedCandidateNodes;
        proposedPowers = base.proposedAssociationPowers;
    otherwise
        error('uav:incm:InvalidTemporalCovfitMethod','未知的时序CovFit方法。');
end
diagnostics.proposedCandidateNodes = proposed;
pairs = associate(priorState, proposed, proposedPowers, settings);
if isempty(pairs)
    return;
end
desired = base.desiredCandidateNodes;
current = proposed(pairs(:,2),:);
transported = priorState.powers(pairs(:,1));
nearDesired = distances(current, desired) <= settings.mergeDistance;
confirmedOutside = false(size(current,1),1);
if resolveConflict
    if isfield(priorState,'confirmedOutside')
        confirmedOutside = priorState.confirmedOutside(pairs(:,1));
    else
        % 兼容旧状态，只能从其可核实的区外位置建立确认。
        confirmedOutside = uav.spectrum.admissibleDirections(priorState.nodes(pairs(:,1),:),cfg);
    end
    % 第一行是固定名义节点，第二行是本期前端临时提名。只有已在区外确认过的
    % 一条轨迹能接管临时峰，名义节点及多候选冲突仍采用保守规则。
    if size(desired,1) > 1 && nnz(nearDesired(:,2)) == 1 && ...
            any(confirmedOutside & nearDesired(:,2) & ~nearDesired(:,1))
        diagnostics.reclassifiedDesiredPeak = desired(2,:);
        desired(2,:) = [];
        nearDesired(:,2) = [];
    end
end
duplicateDesired = any(nearDesired,2);
current(duplicateDesired,:) = [];
confirmedOutside(duplicateDesired) = [];
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
nextState = makeState(nodes, powers, cfg, resolveConflict, ...
    [confirmedOutside;false(size(fresh,1),1)]);
diagnostics.candidateManifold = A;
diagnostics.desiredCandidateNodes = desired;
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

function state = makeState(nodes, powers, cfg, resolveConflict, inherited)
% 只传递本周期拟合可靠的主干扰。区外确认可随匹配轨迹在区内延续。
reference = max([powers;0]);
keep = powers > cfg.incm.temporalCovfit.minimumTrackPowerRatio * reference;
state = struct('nodes', nodes(keep,:), 'powers', powers(keep));
if resolveConflict
    confirmed = uav.spectrum.admissibleDirections(state.nodes,cfg);
    if nargin >= 5
        confirmed = confirmed | inherited(keep);
    end
    state.confirmedOutside = confirmed;
end
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
