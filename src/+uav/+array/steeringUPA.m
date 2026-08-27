function a = steeringUPA(xi, Nx, Ny, dOverLambda)
%STEERINGUPA Return UPA steering vectors in vx-kron-vy sensor order.

if ~isnumeric(xi) || size(xi, 1) ~= 2 || any(~isfinite(xi), 'all')
    error('uav:array:InvalidDirectionCosines', 'xi must be finite with two rows.');
end
if ~isscalar(Nx) || Nx <= 0 || Nx ~= floor(Nx) || ...
        ~isscalar(Ny) || Ny <= 0 || Ny ~= floor(Ny)
    error('uav:array:InvalidDimensions', 'Nx and Ny must be positive integers.');
end
if isscalar(dOverLambda)
    dOverLambda = [dOverLambda, dOverLambda];
end
if ~isnumeric(dOverLambda) || numel(dOverLambda) ~= 2 || ...
        any(~isfinite(dOverLambda)) || any(dOverLambda <= 0)
    error('uav:array:InvalidSpacing', 'dOverLambda must contain two positive finite spacings.');
end
dOverLambda = reshape(dOverLambda, 1, []);

vx = exp(1i * 2 * pi * dOverLambda(1) * xi(1, :) .* (0:Nx - 1).');
vy = exp(1i * 2 * pi * dOverLambda(2) * xi(2, :) .* (0:Ny - 1).');
a = zeros(Nx * Ny, size(xi, 2));
for sourceIndex = 1:size(xi, 2)
    a(:, sourceIndex) = kron(vx(:, sourceIndex), vy(:, sourceIndex));
end
end
