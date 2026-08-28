function [x, probability] = finiteCdf(values)
%FINITECDF 计算仅含有限样本的经验累积分布函数坐标。
%   x为升序样本，probability(i)=i/N。空输入返回两个0×1列向量。
%   本函数不创建图窗，也不依赖Statistics Toolbox的ecdf状态。

x = sort(values(isfinite(values)));
x = x(:);
if isempty(x)
    probability = zeros(0, 1);
else
    probability = (1:numel(x)).' / numel(x);
end
end
