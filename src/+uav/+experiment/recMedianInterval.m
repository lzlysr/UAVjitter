function bounds=recMedianInterval(values,confidence)
%RECMEDIANINTERVAL 二项分布次序统计量给出的保守中位数区间，不假设正态。
%   样本不足以形成有限双侧区间时返回NaN，不伪造确定性。
values=sort(values(isfinite(values))); n=numel(values); bounds=[NaN,NaN];
if n==0,return;end
j=0:n; probabilities=exp(gammaln(n+1)-gammaln(j+1)-gammaln(n-j+1)-n*log(2));
tails=cumsum(probabilities); index=find(tails<=(1-confidence)/2,1,'last');
if ~isempty(index) && index<=n, bounds=[values(index),values(n-index+1)];end
end
