function scenario = generate(cfg, stream)
%GENERATE Generate stationary attitude and steering truth histories.

validateStream(stream);
K = cfg.timeline.K;
M = cfg.derived.M;
sourceCount = cfg.derived.J + 1;
F = diag(cfg.attitude.rho(:));
stdRad = deg2rad(cfg.attitude.stdDeg(:));
Qeta = diag(stdRad .^ 2);
Qnu = Qeta - F * Qeta * F.';
if any(diag(Qnu) < -eps(max(1, max(diag(Qeta)))))
    error('uav:scenario:InvalidAttitudeModel', 'Stationary innovation covariance is not PSD.');
end
innovationStd = sqrt(max(diag(Qnu), 0));

attitudeOffsetRad = zeros(3, K);
attitudeOffsetRad(:, 1) = stdRad .* randn(stream, 3, 1);
for timeIndex = 2:K
    attitudeOffsetRad(:, timeIndex) = F * attitudeOffsetRad(:, timeIndex - 1) + ...
        innovationStd .* randn(stream, 3, 1);
end

nominalRad = deg2rad(cfg.attitude.nominalDeg(:));
attitudeRad = nominalRad + attitudeOffsetRad;
rotations = zeros(3, 3, K);
arrayUnitVectors = zeros(3, sourceCount, K);
directionCosines = zeros(2, sourceCount, K);
azElRad = zeros(2, sourceCount, K);
steeringTruth = complex(zeros(M, sourceCount, K));
for timeIndex = 1:K
    R = uav.geometry.rotationZYX(attitudeRad(:, timeIndex));
    uArray = uav.geometry.worldToArray(R, cfg.derived.worldUnitVectors);
    [xi, directions] = uav.geometry.directionCoordinates(uArray);

    rotations(:, :, timeIndex) = R;
    arrayUnitVectors(:, :, timeIndex) = uArray;
    directionCosines(:, :, timeIndex) = xi;
    azElRad(:, :, timeIndex) = directions;
    steeringTruth(:, :, timeIndex) = uav.array.steeringUPA(xi, cfg.array.Nx, ...
        cfg.array.Ny, cfg.array.spacingOverLambda);
end

scenario = struct( ...
    'F', F, ...
    'Qeta', Qeta, ...
    'Qnu', Qnu, ...
    'attitudeOffsetRad', attitudeOffsetRad, ...
    'attitudeRad', attitudeRad, ...
    'rotations', rotations, ...
    'arrayUnitVectors', arrayUnitVectors, ...
    'directionCosines', directionCosines, ...
    'azElRad', azElRad, ...
    'steeringTruth', steeringTruth);
end

function validateStream(stream)
if ~isa(stream, 'RandStream')
    error('uav:scenario:MissingStream', 'A RandStream input is required.');
end
end
