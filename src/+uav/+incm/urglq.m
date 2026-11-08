function [RinHat,diagnostics]=urglq(obs,cfg)
%URGLQ 仅复现URGLQ局部积分/高斯求积，不使用投影或期望导向校正。
%   论文DOI 10.1109/TAES.2023.3263386，式(19)、(27)、(28)。
%   本对照用原SCM替代论文投影后的准协方差；UPA是明示的二维扩展。
%   不做谱底扣除、功率拟合或结构收缩，也不在积分后重复加入噪声I。
front=uav.spectrum.interferenceCandidates(obs,cfg);
[nodes,weights,quadrature]=uav.spectrum.urglqQuadrature(front.nodes,cfg);
A=uav.array.steeringUPA(nodes.',cfg.array.Nx,cfg.array.Ny,cfg.array.spacingOverLambda);
[powers,powerDiagnostics]=uav.spectrum.caponAtDirections(obs.Rx,A,cfg.incm.urglq.powerLoadingRatio);
if isempty(weights)
    % 没有候选/所有节点被排除时采用明确的噪声阵回退，避免零矩阵求权值。
    RinHat=front.noisePower*eye(cfg.derived.M);fallback="noise-only-empty-domain";
else
    RinHat=(A.*(weights.*powers).')*A';fallback="none";
end
RinHat=(RinHat+RinHat')/2;
diagnostics=struct('method',"incm_urglq",'interfererCandidateNodes',front.nodes, ...
    'nodes',nodes,'weights',weights,'caponPowers',powers,'quadrature',quadrature, ...
    'candidateGridSize',front.candidateGridSize,'candidateLoading',front.candidateLoading, ...
    'peakSelection',front.selection,'powerSpectrum',powerDiagnostics, ...
    'noisePower',front.noisePower,'noiseSource',front.noiseSource,'fallback',fallback, ...
    'desiredProjectionApplied',false,'desiredSteeringCorrectionApplied',false, ...
    'completePaperReproduction',false,'isUpaExtension',cfg.derived.arrayType=="UPA", ...
    'powerFittingApplied',false,'refinementApplied',false,'usedTruth',false);
end
