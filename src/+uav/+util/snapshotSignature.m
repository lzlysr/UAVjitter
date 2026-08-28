function signature = snapshotSignature(X)
%SNAPSHOTSIGNATURE 计算快拍矩阵的确定性身份摘要。
%   摘要依次为实部和、虚部和、Frobenius范数和总能量。它只用于测试
%   “六种方法是否收到同一个X”，不是哈希函数，也不是论文性能指标。

if ~isnumeric(X) || ~ismatrix(X) || isempty(X) || any(~isfinite(X), 'all')
    error('uav:util:InvalidSnapshot', 'X must be a finite nonempty matrix.');
end
signature = [sum(real(X), 'all'), sum(imag(X), 'all'), ...
    norm(X, 'fro'), real(sum(abs(X).^2, 'all'))];
end
