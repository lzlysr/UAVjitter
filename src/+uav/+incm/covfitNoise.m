function [noisePower, noiseSource, noiseDiagnostics] = covfitNoise(Rx,snapshotCount,cfg)
%COVFITNOISE 为Capon和纯子空间CovFit使用相同的小特征值估噪规则。
M=cfg.derived.M;
if snapshotCount<M
    noisePower=cfg.signal.noisePower;
    noiseSource="configured-underdetermined";
    noiseDiagnostics=struct('eigenvalues',zeros(0,1),'lowerEigenvalueCount',0);
else
    % 与原细化版一致：较小特征值的中位数估噪，不改成真值噪声。
    eigenvalues=sort(real(eig((Rx+Rx')/2)),'ascend');
    lowerCount=max(1,ceil(M*cfg.spectrum.noiseEigenvalueFraction));
    rawEstimate=median(eigenvalues(1:lowerCount));
    scale=max(real(trace(Rx))/M,eps);
    noisePower=max(rawEstimate,eps(scale));noiseSource="lower-eigenvalues";
    noiseDiagnostics=struct('eigenvalues',eigenvalues,'lowerEigenvalueCount',lowerCount);
end
end
