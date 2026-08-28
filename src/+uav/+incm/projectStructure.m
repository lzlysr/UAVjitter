function [structured, diagnostics] = projectStructure(R, cfg, applyEnhancement)
%PROJECTSTRUCTURE Hermitianize and optionally floor/shrink a covariance.

if nargin < 3
    applyEnhancement = true;
end
validateInputs(R, cfg, applyEnhancement);
hermitian = (R + R') / 2;
M = size(R, 1);
if ~applyEnhancement
    structured = hermitian;
    diagnostics = struct('applied', false, 'eigenvalueFloor', NaN, ...
        'shrinkage', 0, 'minimumEigenvalue', min(real(eig(hermitian))));
    return;
end

[vectors, values] = eig(hermitian, 'vector');
values = real(values);
covarianceScale = max(real(trace(hermitian)) / M, eps);
eigenvalueFloor = cfg.incm.eigenvalueFloorRatio * covarianceScale;
values = max(values, eigenvalueFloor);
structured = vectors * diag(values) * vectors';
structured = (structured + structured') / 2;
isotropicCovariance = real(trace(structured)) / M * eye(M);
structured = (1 - cfg.incm.shrinkage) * structured + ...
    cfg.incm.shrinkage * isotropicCovariance;
structured = (structured + structured') / 2;
diagnostics = struct( ...
    'applied', true, ...
    'eigenvalueFloor', eigenvalueFloor, ...
    'shrinkage', cfg.incm.shrinkage, ...
    'minimumEigenvalue', min(real(eig(structured))));
end

function validateInputs(R, cfg, applyEnhancement)
if ~ismatrix(R) || size(R, 1) ~= size(R, 2) || any(~isfinite(R), 'all')
    error('uav:incm:InvalidCovariance', 'R must be a finite square covariance.');
end
if ~isstruct(cfg) || ~isfield(cfg, 'incm') || ...
        ~isfield(cfg.incm, 'eigenvalueFloorRatio') || ~isfield(cfg.incm, 'shrinkage')
    error('uav:incm:InvalidConfiguration', 'cfg is missing structure fields.');
end
if ~islogical(applyEnhancement) || ~isscalar(applyEnhancement)
    error('uav:incm:InvalidStructureSwitch', ...
        'applyEnhancement must be a logical scalar.');
end
end
