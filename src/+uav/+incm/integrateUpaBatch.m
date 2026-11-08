function [Rhat,diagnostics] = integrateUpaBatch(Rx,cfg,gridSize,rule,settings)
%INTEGRATEUPABATCH 分块计算多组SCM的同一个二维原始Rec积分，不改变估计公式。
%   Rx为M×M×B；端点规则严格对应当前生产网格及h^2权重。
%   中点规则在同一[-1,1]^2区域错开半个单元，积分测度仍为dxi_x dxi_y。
%   只按可见圆盘及固定期望排除区筛选节点，不读取真实方向。
%   利用UPA重复阵元差，将a^H*Rx^{-1}*a化成有限三角多项式；
%   再计算同一批三角矩。此为代数等价加速，不加载、不插值、不拟合功率。
if cfg.array.Nx<=1 || cfg.array.Ny<=1
    error('uav:incm:UpaOnly','本网格诊断只支持二维UPA。');
end
M=cfg.derived.M; B=size(Rx,3);
if size(Rx,1)~=M || size(Rx,2)~=M || any(~isfinite(Rx),'all')
    error('uav:incm:InvalidBatch','协方差尺寸不匹配或包含非有限值。');
end
if rule=="endpoint"
    h=2/(gridSize-1); axisValues=linspace(-1,1,gridSize);
elseif rule=="midpoint"
    h=2/gridSize; axisValues=-1+h/2+(0:gridSize-1)*h;
else
    error('uav:incm:InvalidQuadrature','规则必须为endpoint或midpoint。');
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
moments=zeros(1+2*Q,B); includedCount=0;
for first=1:settings.nodeChunkSize:gridSize^2
    indices=first:min(first+settings.nodeChunkSize-1,gridSize^2);
    % MATLAB meshgrid列优先展开：y先变化。
    nodes=[axisValues(floor((indices-1)/gridSize)+1).', ...
        axisValues(mod(indices-1,gridSize)+1).'];
    keep=uav.spectrum.admissibleDirections(nodes,cfg); nodes=nodes(keep,:);
    if isempty(nodes),continue;end
    phase=2*pi*(lags.*cfg.array.spacingOverLambda)*nodes.';
    basis=[ones(1,size(nodes,1));cos(phase);sin(phase)];
    includedCount=includedCount+size(nodes,1);
    for start=1:settings.matrixBatchSize:B
        selection=start:min(start+settings.matrixBatchSize-1,B);
        denominator=basis.'*coefficients(:,selection);
        if any(denominator<=0,'all') || any(~isfinite(denominator),'all')
            error('uav:incm:InvalidCaponDenominator','批量Capon分母必须为正且有限。');
        end
        moments(:,selection)=moments(:,selection)+basis*(h^2./denominator);
    end
end
Rhat=complex(zeros(M,M,B));
for b=1:B
    R=moments(1,b)*eye(M);
    for q=1:Q
        value=moments(1+q,b)+1i*moments(1+Q+q,b);
        positive=dx==lags(q,1) & dy==lags(q,2);
        negative=dx==-lags(q,1) & dy==-lags(q,2);
        R(positive)=value;R(negative)=conj(value);
    end
    Rhat(:,:,b)=(R+R')/2;
end
diagnostics=struct('gridSize',gridSize,'rule',rule,'step',h, ...
    'includedNodeCount',includedCount,'includedArea',includedCount*h^2, ...
    'measure',"dxi_x_dxi_y",'loading',0,'usedTruth',false, ...
    'algebraicBasisCount',1+2*Q);
end
