function obs = generateGuLeshem2012Trial(cfg, snrDb, snapshotCount, stream)
%GENERATEGULESHEM2012TRIAL 生成Gu–Leshem 2012例1的一次块快拍。
%   每次蒙特卡洛独立抽取期望与干扰的实际DOA；同一块的全部快拍保持
%   该组DOA不变。期望、干扰和噪声均采用独立零均值圆对称复高斯过程。

validateInputs(cfg, snrDb, snapshotCount, stream);
paper = cfg.literature.guLeshem2012;
M = paper.sensorCount;
nominalAnglesDeg = [paper.nominalDesiredDeg, paper.nominalInterfererDeg];
offsetsDeg = paper.lookMismatchHalfWidthDeg * (2 * rand(stream, 1, 3) - 1);
actualAnglesDeg = nominalAnglesDeg + offsetsDeg;
trueASVs = uav.array.steeringULA(deg2rad(actualAnglesDeg), M, ...
    paper.spacingOverLambda);
nominalSteering = uav.array.steeringULA( ...
    deg2rad(paper.nominalDesiredDeg), M, paper.spacingOverLambda);

[scanAnglesDeg, scanWeights] = complementQuadrature(paper);
scanManifold = uav.array.steeringULA(deg2rad(scanAnglesDeg), M, ...
    paper.spacingOverLambda);

desiredPower = paper.noisePower * 10^(snrDb / 10);
interfererPowers = paper.noisePower * 10.^(paper.interfererInrDb / 10);
desiredSymbols = complexGaussian(stream, 1, snapshotCount);
interfererSymbols = complexGaussian(stream, 2, snapshotCount);
noise = sqrt(paper.noisePower) * complexGaussian(stream, M, snapshotCount);
X = sqrt(desiredPower) * trueASVs(:, 1) * desiredSymbols + ...
    trueASVs(:, 2:end) * (sqrt(interfererPowers(:)) .* interfererSymbols) + noise;
Rx = X * X' / snapshotCount;
trueIncm = paper.noisePower * eye(M) + trueASVs(:, 2:end) * ...
    diag(interfererPowers) * trueASVs(:, 2:end)';

truth = struct( ...
    'trueASVs', trueASVs, ...
    'trueIncm', trueIncm, ...
    'desiredPower', desiredPower, ...
    'interfererPowers', interfererPowers);
obs = struct( ...
    'X', X, ...
    'Rx', (Rx + Rx') / 2, ...
    'nominalSteering', nominalSteering, ...
    'nominalAnglesDeg', nominalAnglesDeg, ...
    'actualAnglesDeg', actualAnglesDeg, ...
    'scanAnglesDeg', scanAnglesDeg, ...
    'scanManifold', scanManifold, ...
    'scanWeights', scanWeights, ...
    'snrDb', snrDb, ...
    'snapshotCount', snapshotCount, ...
    'truth', truth);
end

function [anglesDeg, weights] = complementQuadrature(paper)
%COMPLEMENTQUADRATURE 对[-90,0]和[10,90]分别构造梯形求积。
step = paper.angleGridStepDeg;
sector = paper.desiredSectorDeg;
left = linspace(-90, sector(1), round((sector(1) + 90) / step) + 1);
right = linspace(sector(2), 90, round((90 - sector(2)) / step) + 1);
leftWeights = trapezoidWeights(numel(left), deg2rad(step));
rightWeights = trapezoidWeights(numel(right), deg2rad(step));
anglesDeg = [left, right];
weights = [leftWeights, rightWeights];
end

function weights = trapezoidWeights(count, stepRad)
weights = stepRad * ones(1, count);
weights([1, end]) = stepRad / 2;
end

function values = complexGaussian(stream, rows, columns)
%COMPLEXGAUSSIAN 生成单位方差CN(0,1)样本。
values = (randn(stream, rows, columns) + ...
    1i * randn(stream, rows, columns)) / sqrt(2);
end

function validateInputs(cfg, snrDb, snapshotCount, stream)
if ~isstruct(cfg) || ~isfield(cfg, 'literature') || ...
        ~isfield(cfg.literature, 'guLeshem2012')
    error('uav:scenario:MissingGuLeshemConfiguration', ...
        'cfg.literature.guLeshem2012 is required.');
end
if ~isscalar(snrDb) || ~isfinite(snrDb)
    error('uav:scenario:InvalidGuLeshemSnr', 'snrDb must be finite.');
end
if ~isscalar(snapshotCount) || snapshotCount < 1 || ...
        snapshotCount ~= floor(snapshotCount)
    error('uav:scenario:InvalidGuLeshemSnapshots', ...
        'snapshotCount must be a positive integer.');
end
if ~isa(stream, 'RandStream')
    error('uav:scenario:InvalidRandomStream', ...
        'stream must be an explicit RandStream.');
end
end
