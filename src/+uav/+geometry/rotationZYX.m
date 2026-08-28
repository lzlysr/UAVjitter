function R = rotationZYX(etaRad)
%ROTATIONZYX 根据滚转、俯仰和偏航角生成三维旋转矩阵。
%   etaRad = [roll; pitch; yaw]，单位为弧度。
%   输出R满足 v_world = R*v_array，旋转顺序固定为
%   R = Rz(yaw)*Ry(pitch)*Rx(roll)。改变乘法顺序会改变物理含义。

etaRad = etaRad(:);
if numel(etaRad) ~= 3 || any(~isfinite(etaRad))
    error('uav:geometry:InvalidAttitude', 'etaRad must contain three finite radians.');
end

roll = etaRad(1);
pitch = etaRad(2);
yaw = etaRad(3);
% 分别绕阵列坐标系x、y、z轴建立右手系基本旋转。
Rx = [1, 0, 0; 0, cos(roll), -sin(roll); 0, sin(roll), cos(roll)];
Ry = [cos(pitch), 0, sin(pitch); 0, 1, 0; -sin(pitch), 0, cos(pitch)];
Rz = [cos(yaw), -sin(yaw), 0; sin(yaw), cos(yaw), 0; 0, 0, 1];
R = Rz * Ry * Rx;
end
