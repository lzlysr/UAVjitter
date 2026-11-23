function powers = candidateAssociationPowers(A,Rx,noisePower)
%CANDIDATEASSOCIATIONPOWERS 估计各方向候选的功率，仅用于时序关联。
%   两种前端共用双边最小二乘。最终干扰功率仍由协方差拟合给出。
powers = zeros(0,1);
if isempty(A) || rcond(A'*A) <= sqrt(eps)
    return;
end
inverseA = A \ eye(size(Rx,1));
sourceCovariance = inverseA*(Rx-noisePower*eye(size(Rx,1)))*inverseA';
powers = max(real(diag(sourceCovariance)),0);
end
