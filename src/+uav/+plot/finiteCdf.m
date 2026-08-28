function [x, probability] = finiteCdf(values)
%FINITECDF Return a finite-only empirical CDF without toolbox state.

x = sort(values(isfinite(values)));
x = x(:);
if isempty(x)
    probability = zeros(0, 1);
else
    probability = (1:numel(x)).' / numel(x);
end
end
