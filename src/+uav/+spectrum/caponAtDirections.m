function [powers,diagnostics]=caponAtDirections(Rx,A,loadingRatio)
%CAPONATDIRECTIONS 在指定流形列上读取Capon值，绝不引入积分面积或拟合。
%   加载比例0代表原始SCM公式；秩亏时明确报不适用，不静默使用伪逆。
%   正加载必须由config显式指定，加载量为比例乘平均阵元功率。
M=size(Rx,1);Rx=(Rx+Rx')/2;
loading=loadingRatio*max(real(trace(Rx))/M,eps);
R=Rx+loading*eye(M);
if rcond(R)<=eps(class(R))
    error('uav:incm:SingularSampleCovariance','未加载Capon需要可逆SCM；请增加快拍或显式配置加载。');
end
solved=R\A;denominator=real(sum(conj(A).*solved,1));
if any(~isfinite(denominator)) || any(denominator<=0)
    error('uav:spectrum:InvalidCaponPower','指定方向的Capon分母必须为正且有限。');
end
powers=reshape(1./denominator,[],1);
diagnostics=struct('loadingRatio',loadingRatio,'loading',loading, ...
    'reciprocalCondition',rcond(R),'usedPseudoInverse',false);
end
