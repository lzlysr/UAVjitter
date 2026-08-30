function grid = upaIntegrationGrid(cfg, gridSize, measureName)
%UPAINTEGRATIONGRID 生成UPA二维诊断所需的节点、流形和求积权重。
%   节点始终位于名义姿态下的下半球可见域。includedMask仅标记
%   排除期望区后可用于干扰积分的节点；本函数不读取任何真实量。

validateInputs(cfg, gridSize, measureName);
measureName = string(measureName);
switch measureName
    case "direction_cosine"
        [nodes, weights] = directionCosineGrid(gridSize);
        % 主积分保持方向余弦面积元；公开立体角诊断改由独立等立体角网格计算。
        [coverageNodes, coverageWeights] = azimuthZGrid(gridSize);
    case "solid_angle"
        [nodes, weights] = azimuthZGrid(gridSize);
        coverageNodes = nodes;
        coverageWeights = weights;
    case "az_el_solid_angle"
        [nodes, weights] = azElGrid(gridSize);
        coverageNodes = nodes;
        coverageWeights = weights;
end

manifold = uav.array.steeringUPA(nodes.', cfg.array.Nx, cfg.array.Ny, ...
    cfg.array.spacingOverLambda);
admissibleMask = uav.spectrum.admissibleDirections(nodes, cfg);
desiredMask = ~admissibleMask;
includedMask = ~desiredMask;
coverageMask = uav.spectrum.admissibleDirections(coverageNodes, cfg);
grid = struct( ...
    'nodes', nodes, ...
    'manifold', manifold, ...
    'weights', reshape(weights, 1, []), ...
    'representedSolidAngle', sum(coverageWeights(coverageMask)), ...
    'measureName', measureName, ...
    'desiredMask', desiredMask, ...
    'includedMask', includedMask);
end

function validateInputs(cfg, gridSize, measureName)
%VALIDATEINPUTS 校验已完成配置派生和离散求积的必要输入。
if ~isstruct(cfg) || ~isfield(cfg, 'array') || ~isfield(cfg, 'spectrum') || ...
        ~isfield(cfg, 'derived') || ~isfield(cfg.derived, 'nominalArrayUnitVectors')
    error('uav:spectrum:InvalidConfiguration', ...
        'A validated configuration with nominal array directions is required.');
end
if ~isnumeric(gridSize) || ~isscalar(gridSize) || ~isfinite(gridSize) || ...
        gridSize < 3 || gridSize ~= floor(gridSize) || mod(gridSize, 2) == 0
    error('uav:spectrum:InvalidGridSize', ...
        'gridSize must be an odd integer greater than or equal to three.');
end
allowedMeasures = ["direction_cosine", "solid_angle", "az_el_solid_angle"];
if ~(ischar(measureName) || (isstring(measureName) && isscalar(measureName))) || ...
        ~ismember(string(measureName), allowedMeasures)
    error('uav:spectrum:InvalidIntegrationMeasure', ...
        'measureName must name a supported UPA integration measure.');
end
end

function [nodes, weights] = directionCosineGrid(gridSize)
%DIRECTIONCOSINEGRID 在方向余弦单元中心构造可见圆盘上的节点。
step = 2 / gridSize;
axisValues = -1 + step / 2 + step * (0:gridSize - 1);
[xiX, xiY] = meshgrid(axisValues, axisValues);
visibleMask = xiX.^2 + xiY.^2 <= 1;
nodes = [xiX(visibleMask), xiY(visibleMask)];
weights = step^2 * ones(size(nodes, 1), 1);
end

function [nodes, weights] = azimuthZGrid(gridSize)
%AZIMUTHZGRID 用方位角和z方向余弦的等立体角单元覆盖下半球。
%   z取[-1,0]单元中心，dOmega=daz*dz，不经过可见圆盘边界的
%   奇异Jacobian；因此所有网格大小的总权重严格等于2*pi。
daz = 2 * pi / gridSize;
dz = 1 / gridSize;
azimuth = -pi + daz / 2 + daz * (0:gridSize - 1);
z = -1 + dz / 2 + dz * (0:gridSize - 1);
[azimuthGrid, zGrid] = meshgrid(azimuth, z);
radial = sqrt(1 - zGrid(:).^2);
nodes = [radial .* cos(azimuthGrid(:)), ...
    radial .* sin(azimuthGrid(:))];
weights = daz * dz * ones(gridSize^2, 1);
end

function [nodes, weights] = azElGrid(gridSize)
%AZELGRID 在方位角和下半球高度角单元中心构造立体角求积节点。
daz = 2 * pi / gridSize;
del = pi / (2 * gridSize);
azimuth = -pi + daz / 2 + daz * (0:gridSize - 1);
elevation = -pi / 2 + del / 2 + del * (0:gridSize - 1);
[azimuthGrid, elevationGrid] = meshgrid(azimuth, elevation);
nodes = [cos(elevationGrid(:)) .* cos(azimuthGrid(:)), ...
    cos(elevationGrid(:)) .* sin(azimuthGrid(:))];
weights = cos(elevationGrid(:)) * daz * del;
end
