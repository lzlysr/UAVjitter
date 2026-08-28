function c = temporalProx(p, g, lambdaS, lambdaT)
%TEMPORALPROX Exact nonnegative prox for two absolute-value penalties.
%
% Minimizes (p-c)^2 + lambdaS*c + lambdaT*abs(c-g) independently for
% every element, with p >= 0, g >= 0, and c >= 0.

validateInputs(p, g, lambdaS, lambdaT);
combinedInput = p + g + lambdaS + lambdaT;
outputSize = size(combinedInput);
p = p + zeros(outputSize);
g = g + zeros(outputSize);
lambdaS = lambdaS + zeros(outputSize);
lambdaT = lambdaT + zeros(outputSize);

pVector = p(:);
gVector = g(:);
lambdaSVector = lambdaS(:);
lambdaTVector = lambdaT(:);
betweenCandidate = min(max(pVector - (lambdaSVector - lambdaTVector) / 2, 0), ...
    gVector);
aboveCandidate = max(pVector - (lambdaSVector + lambdaTVector) / 2, gVector);
candidates = [zeros(numel(pVector), 1), gVector, betweenCandidate, aboveCandidate];
objective = (pVector - candidates).^2 + lambdaSVector .* candidates + ...
    lambdaTVector .* abs(candidates - gVector);
[~, selected] = min(objective, [], 2);
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
