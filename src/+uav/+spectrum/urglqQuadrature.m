function [nodes,weights,diagnostics]=urglqQuadrature(candidates,cfg)
%URGLQQUADRATURE 在干扰局部域构造高斯节点；节点不是独立物理信号源。
%   ULA：等效角theta=asin(u_x)，局部区间合并后积分dtheta。
%   UPA：负半球上的方位/高度矩形，跨正负pi时拆分；测度为dOmega。
%   重叠矩形按所有边界分割为不重叠单元，避免多峰邻域重复计入功率。
%   固定期望排除区仍在节点处执行，边界裁剪下低阶求积只是近似。
s=cfg.incm.urglq;[z,w]=gaussLegendre(s.quadratureOrder);
nodes=zeros(0,2);weights=zeros(0,1);domainCount=0;
if cfg.derived.arrayType=="ULA"
    theta=asin(max(-1,min(1,candidates(:,1))));h=deg2rad(s.ulaHalfWidthDeg);
    intervals=sortrows([max(-pi/2,theta-h),min(pi/2,theta+h)],1);
    merged=zeros(0,2);
    for k=1:size(intervals,1)
        if isempty(merged) || intervals(k,1)>merged(end,2)
            merged(end+1,:)=intervals(k,:); %#ok<AGROW>
        else
            merged(end,2)=max(merged(end,2),intervals(k,2));
        end
    end
    for k=1:size(merged,1)
        center=mean(merged(k,:));half=diff(merged(k,:))/2;
        angle=center+half*z;nodes=[nodes;sin(angle),zeros(size(angle))]; %#ok<AGROW>
        weights=[weights;half*w];domainCount=domainCount+1; %#ok<AGROW>
    end
    measure="dtheta";
else
    half=deg2rad(s.upaAzElHalfWidthDeg);rectangles=zeros(0,4);
    azimuth=atan2(candidates(:,2),candidates(:,1));
    elevation=asin(-sqrt(max(0,1-sum(candidates.^2,2))));
    for k=1:size(candidates,1)
        % 周期方向先平移后裁剪，保证180/-180度附近的积分域连续。
        for shift=[-2*pi,0,2*pi]
            left=max(-pi,azimuth(k)-half(1)+shift);right=min(pi,azimuth(k)+half(1)+shift);
            if left<right
                rectangles(end+1,:)=[left,right,max(-pi/2,elevation(k)-half(2)),min(0,elevation(k)+half(2))]; %#ok<AGROW>
            end
        end
    end
    x=unique(reshape(rectangles(:,1:2),[],1));y=unique(reshape(rectangles(:,3:4),[],1));
    [Zx,Zy]=meshgrid(z,z);[Wx,Wy]=meshgrid(w,w);
    for ix=1:numel(x)-1
        for iy=1:numel(y)-1
            center=[mean(x(ix:ix+1)),mean(y(iy:iy+1))];
            inside=rectangles(:,1)<center(1) & rectangles(:,2)>center(1) & ...
                rectangles(:,3)<center(2) & rectangles(:,4)>center(2);
            if ~any(inside),continue;end
            radius=[diff(x(ix:ix+1)),diff(y(iy:iy+1))]/2;
            az=center(1)+radius(1)*Zx(:);el=center(2)+radius(2)*Zy(:);
            nodes=[nodes;cos(el).*cos(az),cos(el).*sin(az)]; %#ok<AGROW>
            % dOmega=cos(el)daz del。3阶时单矩形为3×3节点。
            weights=[weights;prod(radius)*Wx(:).*Wy(:).*cos(el)]; %#ok<AGROW>
            domainCount=domainCount+1;
        end
    end
    measure="solid_angle";
end
beforeCount=size(nodes,1);beforeArea=sum(weights);
keep=uav.spectrum.admissibleDirections(nodes,cfg);nodes=nodes(keep,:);weights=weights(keep);
diagnostics=struct('order',s.quadratureOrder,'measure',measure,'domainCellCount',domainCount, ...
    'nodeCountBeforeDesiredExclusion',beforeCount,'nodeCount',size(nodes,1), ...
    'areaBeforeDesiredExclusion',beforeArea,'representedArea',sum(weights), ...
    'overlapPolicy',"nonoverlapping-union",'usedTruth',false);
end

function [nodes,weights]=gaussLegendre(order)
% Golub–Welsch构造[-1,1]节点及权重，3阶即±sqrt(3/5)、0和5/9、8/9、5/9。
k=(1:order-1).';off=k./sqrt(4*k.^2-1);
J=diag(off,1)+diag(off,-1);
if order==1,J=0;end
[V,D]=eig(J);[nodes,index]=sort(diag(D));weights=2*(V(1,index).^2).';
end
