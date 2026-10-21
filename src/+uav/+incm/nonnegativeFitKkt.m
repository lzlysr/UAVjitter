function residual = nonnegativeFitKkt(Rx,A,noise,p,lambda)
%NONNEGATIVEFITKKT 独立检查协方差NNLS的非负、对偶可行性和互补条件。
%   返回相对KKT违反程度，不修正返回功率，也不改变生产拟合求解器。
H=abs(A'*A).^2+lambda*eye(size(A,2));
b=real(sum(conj(A).*((Rx-noise*eye(size(Rx,1)))*A),1)).';
g=H*p-b;scale=max([norm(H,inf)*norm(p,inf),norm(b,inf),eps]);
residual=max([max([0;-g])/scale, ...
    norm(p.*g,inf)/(max(norm(p,inf),eps)*scale), ...
    max([0;-p])/max(norm(p,inf),eps)]);
if ~isfinite(residual),residual=Inf;end
end
