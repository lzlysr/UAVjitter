function [Rhat,diagnostics] = integrateUpaBatch(Rx,cfg,gridSize,rule,settings,measures)
%INTEGRATEUPABATCH 分块计算多组SCM的同一个二维原始Rec积分，不改变估计公式。
%   Rx为M×M×B；端点规则严格对应当前生产网格及h^2权重。
%   中点规则在同一[-1,1]^2区域错开半个单元，积分测度仍为dxi_x dxi_y。
%   只按可见圆盘及固定期望排除区筛选节点，不读取真实方向。
%   利用UPA重复阵元差，将a^H*Rx^{-1}*a化成有限三角多项式；
%   再计算同一批三角矩。此为代数等价加速，不加载、不插值、不拟合功率。
%   研究分支新增surface_*规则：以方位phi和负z坐标等距求积。
%   dOmega=dphi dz，dxi_x dxi_y=|z|dphi dz；两测度共用所有节点和Capon值。
%   不在圆盘边缘直接除以零，也不对Jacobian人为裁剪。
%   多测度时Rhat第四维对应measures，diagnostics也按测度排列。
%   azel_*规则直接均匀采样方位/高度角；angle_product表示dphi dalpha，
%   仅作为另一种估计器的诊断，不能把它与dOmega无条件视为同一功率密度。
if nargin<6,measures="dxi_x_dxi_y";end
measures=reshape(string(measures),1,[]);measureCount=numel(measures);
if ~all(ismember(measures,["dxi_x_dxi_y","solid_angle","angle_product"]))
    error('uav:incm:InvalidMeasure','仅支持方向余弦面积、立体角及角度乘积测度。');
end
if cfg.array.Nx<=1 || cfg.array.Ny<=1
    error('uav:incm:UpaOnly','本网格诊断只支持二维UPA。');
end
M=cfg.derived.M; B=size(Rx,3);
if size(Rx,1)~=M || size(Rx,2)~=M || any(~isfinite(Rx),'all')
    error('uav:incm:InvalidBatch','协方差尺寸不匹配或包含非有限值。');
end
surfaceRule=ismember(rule,["surface_midpoint","surface_shifted"]);
azelRule=ismember(rule,["azel_midpoint","azel_shifted"]);
if any(measures=="angle_product") && ~azelRule
    error('uav:incm:InvalidAngleRule','角度乘积测度必须使用azel_*规则。');
end
if rule=="endpoint"
    h=2/(gridSize-1); axisValues=linspace(-1,1,gridSize);
elseif rule=="midpoint"
    h=2/gridSize; axisValues=-1+h/2+(0:gridSize-1)*h;
elseif surfaceRule || azelRule
    if azelRule,h=[2*pi/gridSize,pi/(2*gridSize)];else,h=[2*pi/gridSize,1/gridSize];end
    phi=-pi+h(1)/2+(0:gridSize-1)*h(1);
    if ismember(rule,["surface_shifted","azel_shifted"]),phi=phi-h(1)/2;end
    if azelRule
        elevation=-pi/2+h(2)/2+(0:gridSize-1)*h(2);zAxis=sin(elevation);
    else
        zAxis=-1+h(2)/2+(0:gridSize-1)*h(2);
    end
else
    error('uav:incm:InvalidQuadrature','不支持的求积规则。');
end
if ~surfaceRule && ~azelRule && any(measures~="dxi_x_dxi_y")
    error('uav:incm:SingularBoundaryRule','立体角积分请使用surface_*规则避免边缘奇异性。');
end
% 阵元顺序与steeringUPA的vx kron vy相同。只保存正半平面的阵元差。
x=kron((0:cfg.array.Nx-1)',ones(cfg.array.Ny,1));
y=repmat((0:cfg.array.Ny-1)',cfg.array.Nx,1);
dx=x-x.';dy=y-y.';
lags=unique([dx(:),dy(:)],'rows');
lags=lags(lags(:,1)>0 | (lags(:,1)==0 & lags(:,2)>0),:);
Q=size(lags,1); coefficients=zeros(1+2*Q,B);
for b=1:B
    covariance=(Rx(:,:,b)+Rx(:,:,b)')/2;
    if rcond(covariance)<=eps
        error('uav:incm:SingularSampleCovariance','原始Rec要求可逆协方差。');
    end
    inverse=covariance\eye(M);
    coefficients(1,b)=real(trace(inverse));
    for q=1:Q
        z=sum(inverse(dx==lags(q,1) & dy==lags(q,2)));
        coefficients(1+q,b)=2*real(z);
        coefficients(1+Q+q,b)=2*imag(z);
    end
end
moments=zeros(1+2*Q,B,measureCount); includedCount=0;includedArea=zeros(1,measureCount);
% 可选本地并行求和只用于大规模诊断；默认仍是可独立运行的MATLAB内核。
useNative=isfield(settings,'nativeAcceleration') && settings.nativeAcceleration;
if useNative
    if exist('upaPopulationQuadratureNative','file')~=3
        error('uav:incm:MissingNativeKernel','请先编译诊断MEX并加入路径，或关闭nativeAcceleration。');
    end
    if surfaceRule || azelRule,firstAxis=phi;secondAxis=zAxis;
    else,firstAxis=axisValues;secondAxis=axisValues;end
    if rule=="endpoint",ruleId=0;elseif rule=="midpoint",ruleId=1;
    elseif surfaceRule,ruleId=2;else,ruleId=3;end
    ids=zeros(1,measureCount);ids(measures=="dxi_x_dxi_y")=1;
    ids(measures=="solid_angle")=2;ids(measures=="angle_product")=3;
    [moments,includedCount,includedArea]=upaPopulationQuadratureNative( ...
        coefficients,lags,firstAxis,secondAxis,ruleId, ...
        [cfg.array.Nx,cfg.array.Ny,cfg.array.spacingOverLambda], ...
        cfg.derived.nominalArrayUnitVectors(:,1), ...
        deg2rad(cfg.spectrum.desiredExclusionHalfWidthDeg),ids,settings.nativeThreadCount);
else
for first=1:settings.nodeChunkSize:gridSize^2
    indices=first:min(first+settings.nodeChunkSize-1,gridSize^2);
    if surfaceRule || azelRule
        % z坐标先变化，方位角是周期坐标；没有重复的正负pi节点。
        az=phi(floor((indices-1)/gridSize)+1).';
        z=zAxis(mod(indices-1,gridSize)+1).';r=sqrt(1-z.^2);
        nodes=[r.*cos(az),r.*sin(az)];
    else
        % MATLAB meshgrid列优先展开：y先变化。
        nodes=[axisValues(floor((indices-1)/gridSize)+1).', ...
            axisValues(mod(indices-1,gridSize)+1).'];
    end
    keep=uav.spectrum.admissibleDirections(nodes,cfg); nodes=nodes(keep,:);
    if isempty(nodes),continue;end
    weights=zeros(size(nodes,1),measureCount);
    for m=1:measureCount
        if azelRule
            % dOmega=cos(alpha)dphi dalpha；du dv=|sin(alpha)|dOmega。
            if measures(m)=="angle_product"
                weights(:,m)=prod(h);
            elseif measures(m)=="solid_angle"
                weights(:,m)=r(keep)*prod(h);
            else
                weights(:,m)=r(keep).*abs(z(keep))*prod(h);
            end
        elseif ~surfaceRule
            weights(:,m)=h^2;
        elseif measures(m)=="solid_angle"
            weights(:,m)=prod(h);
        else
            weights(:,m)=abs(z(keep))*prod(h);
        end
    end
    includedArea=includedArea+sum(weights,1);
    % 用整数阵元差的幂递推代替逐阵元差求sin/cos，仍是同一个有限傅里叶基。
    % 9×9阵列只需两个基本复指数，不改变节点、Capon值或求积权重。
    basis=lagBasis(nodes,lags,cfg);
    includedCount=includedCount+size(nodes,1);
    for start=1:settings.matrixBatchSize:B
        selection=start:min(start+settings.matrixBatchSize-1,B);
        denominator=basis.'*coefficients(:,selection);
        if any(denominator<=0,'all') || any(~isfinite(denominator),'all')
            error('uav:incm:InvalidCaponDenominator','批量Capon分母必须为正且有限。');
        end
        for m=1:measureCount
            moments(:,selection,m)=moments(:,selection,m)+basis*(weights(:,m)./denominator);
        end
    end
end
end
Rhat=complex(zeros(M,M,B,measureCount));
for m=1:measureCount
    for b=1:B
        R=moments(1,b,m)*eye(M);
        for q=1:Q
            value=moments(1+q,b,m)+1i*moments(1+Q+q,b,m);
            positive=dx==lags(q,1) & dy==lags(q,2);
            negative=dx==-lags(q,1) & dy==-lags(q,2);
            R(positive)=value;R(negative)=conj(value);
        end
        Rhat(:,:,b,m)=(R+R')/2;
    end
    diagnostics(m)=struct('gridSize',gridSize,'rule',rule,'step',h, ...
        'includedNodeCount',includedCount,'includedArea',includedArea(m), ...
        'measure',measures(m),'loading',0,'usedTruth',false, ...
        'algebraicBasisCount',1+2*Q,'nativeAcceleration',useNative); %#ok<AGROW>
end
end

function basis=lagBasis(nodes,lags,cfg)
count=size(nodes,1);
px=ones(cfg.array.Nx,count);py=ones(cfg.array.Ny,count);
ex=exp(1i*2*pi*cfg.array.spacingOverLambda(1)*nodes(:,1).');
ey=exp(1i*2*pi*cfg.array.spacingOverLambda(2)*nodes(:,2).');
for k=2:cfg.array.Nx,px(k,:)=px(k-1,:).*ex;end
for k=2:cfg.array.Ny,py(k,:)=py(k-1,:).*ey;end
y=py(abs(lags(:,2))+1,:);
negative=lags(:,2)<0;y(negative,:)=conj(y(negative,:));
values=px(lags(:,1)+1,:).*y;
basis=[ones(1,count);real(values);imag(values)];
end
