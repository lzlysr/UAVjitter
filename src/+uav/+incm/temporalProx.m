function c = temporalProx(p, g, lambdaS, lambdaT)
%TEMPORALPROX 精确求解逐元素非负时序稀疏子问题。
%
%   对每个谱单元独立最小化：
%     (p-c)^2 + lambdaS*c + lambdaT*|c-g|,  c>=0。
%   p为当前候选谱，g为对齐后的上一周期增强谱。非负约束下|c|=c。
%   目标在断点c=0和c=g两侧为分段凸二次函数，因此只需比较有限个
%   解析候选点，不需要CVX或迭代优化器。

validateInputs(p, g, lambdaS, lambdaT);
combinedInput = p + g + lambdaS + lambdaT;
% 利用隐式扩展统一标量正则项和矩阵输入的尺寸。
outputSize = size(combinedInput);
p = p + zeros(outputSize);
g = g + zeros(outputSize);
lambdaS = lambdaS + zeros(outputSize);
lambdaT = lambdaT + zeros(outputSize);

pVector = p(:);
gVector = g(:);
lambdaSVector = lambdaS(:);
lambdaTVector = lambdaT(:);
% 区间0<=c<=g和区间c>=g各有一个驻点，再加两个分段边界0、g。
betweenCandidate = min(max(pVector - (lambdaSVector - lambdaTVector) / 2, 0), ...
    gVector);
aboveCandidate = max(pVector - (lambdaSVector + lambdaTVector) / 2, gVector);
candidates = [zeros(numel(pVector), 1), gVector, betweenCandidate, aboveCandidate];
objective = (pVector - candidates).^2 + lambdaSVector .* candidates + ...
    lambdaTVector .* abs(candidates - gVector);
[~, selected] = min(objective, [], 2);
% 每个谱单元独立选择目标函数最小的候选点。
cVector = candidates(sub2ind(size(candidates), (1:numel(pVector)).', selected));
c = reshape(cVector, outputSize);
end

function validateInputs(p, g, lambdaS, lambdaT)
inputs = {p, g, lambdaS, lambdaT};
for index = 1:numel(inputs)
    value = inputs{index};
    if ~isnumeric(value) || isempty(value) || any(~isfinite(value), 'all') || ...
            any(value < 0, 'all')
        error('uav:incm:InvalidTemporalProxInput', ...
            'Prox inputs must be finite, numeric, nonnegative, and nonempty.');
    end
end
try
    combinedInput = p + g + lambdaS + lambdaT; %#ok<NASGU>
catch
    error('uav:incm:InvalidTemporalProxInput', ...
        'Prox inputs must have compatible sizes.');
end
end
