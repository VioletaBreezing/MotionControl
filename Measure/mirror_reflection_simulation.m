%% 镜面反射与检波面交点仿真（含刚体）
% 场景：刚体（长方体）上有两个任意朝向的平面镜 M1, M2
%       两束独立光线分别入射，各自有检波面，反射后再次穿过检波面

clear; clc; close all;

%% 1. 刚体参数（本体系）
Lx = 2; Ly = 2; Lz = 1;
r1_0 = [-Lx/2; 0; 0];   % M1 安装点（X负向面中心）
r2_0 = [0; -Ly/2; 0];   % M2 安装点（Y负向面中心）

% 镜面法向（本体系，可任意设定）
n1_body = [-1; 0; 0]; n1_body = n1_body / norm(n1_body);
n2_body = [0; -1; 0]; n2_body = n2_body / norm(n2_body);

%% 2. 刚体位姿：小角度旋转 + 平移
% 小角度（弧度）
alpha = deg2rad(0);   % 绕 X 轴（滚转）
beta  = deg2rad(0);   % 绕 Y 轴（俯仰）
gamma = deg2rad(0);   % 绕 Z 轴（偏航）

% 构造旋转矩阵（intrinsic X->Y->Z => R = Rz * Ry * Rx）
Rx = [1, 0, 0;
      0, cos(alpha), -sin(alpha);
      0, sin(alpha),  cos(alpha)];

Ry = [cos(beta), 0, sin(beta);
      0, 1, 0;
      -sin(beta), 0, cos(beta)];

Rz = [cos(gamma), -sin(gamma), 0;
      sin(gamma),  cos(gamma), 0;
      0, 0, 1];

R = Rz * Ry * Rx;   % 总旋转矩阵

t = [0.1; 0.2; 0.3];      % 刚体质心位置

% 镜面在世界系中的参数
P1 = R * r1_0 + t;
N1 = R * n1_body;
P2 = R * r2_0 + t;
N2 = R * n2_body;

%% 3. 入射光线1（对应M1）
a1 = [-5; 0; 0];
u1 = [1; 0; 0]; u1 = u1 / norm(u1);
Sigma1_point = a1;
Sigma1_normal = u1;

%% 4. 入射光线2（对应M2）
a2 = [0; -5; 0];
u2 = [0; 1; 0]; u2 = u2 / norm(u2);
Sigma2_point = a2;
Sigma2_normal = u2;

%% 5. 计算M1光路
denom1 = dot(u1, N1);
if abs(denom1) < 1e-8, error('光线1与M1平行'); end
s1 = dot(P1 - a1, N1) / denom1;
r_inc1 = a1 + s1 * u1;
u1_ref = u1 - 2 * dot(u1, N1) * N1;
denom_u1 = dot(u1_ref, u1);
if abs(denom_u1) > 1e-8
    u_param1 = -s1 / denom_u1;
    r_sigma1 = r_inc1 + u_param1 * u1_ref;
else
    r_sigma1 = NaN(3,1);
end

%% 6. 计算M2光路
denom2 = dot(u2, N2);
if abs(denom2) < 1e-8, error('光线2与M2平行'); end
s2 = dot(P2 - a2, N2) / denom2;
r_inc2 = a2 + s2 * u2;
u2_ref = u2 - 2 * dot(u2, N2) * N2;
denom_u2 = dot(u2_ref, u2);
if abs(denom_u2) > 1e-8
    u_param2 = -s2 / denom_u2;
    r_sigma2 = r_inc2 + u_param2 * u2_ref;
else
    r_sigma2 = NaN(3,1);
end

%% 7. 可视化
figure('Color','w'); hold on; axis equal; grid on;
xlabel('X'); ylabel('Y'); zlabel('Z');
title('双镜反射 + 刚体 + 检波面仿真');

% === 新增：绘制刚体（长方体）===
drawRigidBody(t, R, Lx, Ly, Lz, [0.6 0.6 0.9], 0.2);  % 半透明蓝色刚体

% 绘制镜面（叠加在刚体上，更醒目）
drawPlane(P1, N1, 'r', 1.2, 0.6);   % M1 红色，稍小，半透明
drawPlane(P2, N2, 'b', 1.2, 0.6);   % M2 蓝色

% 光线1
plot3([a1(1), r_inc1(1)], [a1(2), r_inc1(2)], [a1(3), r_inc1(3)], 'k--', 'LineWidth',1.5);
plot3([r_inc1(1), r_sigma1(1)], [r_inc1(2), r_sigma1(2)], [r_inc1(3), r_sigma1(3)], 'k-', 'LineWidth',2);
scatter3(a1(1),a1(2),a1(3), 80, 'ko', 'filled');
scatter3(r_inc1(1),r_inc1(2),r_inc1(3), 80, 'ro', 'filled');
scatter3(r_sigma1(1),r_sigma1(2),r_sigma1(3), 80, 'go', 'filled');

% 光线2
plot3([a2(1), r_inc2(1)], [a2(2), r_inc2(2)], [a2(3), r_inc2(3)], 'm--', 'LineWidth',1.5);
plot3([r_inc2(1), r_sigma2(1)], [r_inc2(2), r_sigma2(2)], [r_inc2(3), r_sigma2(3)], 'm-', 'LineWidth',2);
scatter3(a2(1),a2(2),a2(3), 80, 'mo', 'filled');
scatter3(r_inc2(1),r_inc2(2),r_inc2(3), 80, 'bo', 'filled');
scatter3(r_sigma2(1),r_sigma2(2),r_sigma2(3), 80, 'co', 'filled');

% 检波面
drawPlane(Sigma1_point, Sigma1_normal, [0.8 0.8 0.8], 3, 0.15);
drawPlane(Sigma2_point, Sigma2_normal, [0.9 0.7 0.9], 3, 0.15);

legend({'刚体','M1镜面','M2镜面',...
        '光线1入射','光线1反射',...
        '起点1','M1入射点','交点1',...
        '光线2入射','光线2反射',...
        '起点2','M2入射点','交点2',...
        '检波面1','检波面2'},...
        'Location','bestoutside');

view(35,25);
camlight; lighting gouraud;

%% ================== 辅助函数 ==================

% 绘制刚体（中心在 t，方向由 R 定义，尺寸 Lx,Ly,Lz）
function drawRigidBody(center, R, Lx, Ly, Lz, color, alpha)
    % 8个顶点（本体系）
    dx = Lx/2; dy = Ly/2; dz = Lz/2;
    vertices_body = [
        -dx -dy -dz;
         dx -dy -dz;
         dx  dy -dz;
        -dx  dy -dz;
        -dx -dy  dz;
         dx -dy  dz;
         dx  dy  dz;
        -dx  dy  dz
    ]';
    
    % 变换到世界系：先旋转，再平移
    vertices_world = R * vertices_body + center;
    
    % 定义6个面（每个面4个顶点索引）
    faces = [
        1 2 3 4;  % 底面 (z=-dz)
        5 6 7 8;  % 顶面 (z=+dz)
        1 2 6 5;  % 前面 (y=-dy)
        3 4 8 7;  % 后面 (y=+dy)
        1 4 8 5;  % 左面 (x=-dx)
        2 3 7 6   % 右面 (x=+dx)
    ];
    
    % 绘制每个面
    for i = 1:size(faces,1)
        idx = faces(i,:);
        X = vertices_world(1,idx);
        Y = vertices_world(2,idx);
        Z = vertices_world(3,idx);
        fill3(X,Y,Z, color, 'FaceAlpha',alpha, 'EdgeColor',[0.3 0.3 0.3], 'LineWidth',0.8);
    end
end

% 绘制平面（用于镜面和检波面）
function drawPlane(center, normal, color, size, alpha)
    if nargin < 5, alpha = 0.3; end
    if nargin < 4, size = 1; end
    normal = normal(:) / norm(normal);
    
    % 找两个正交切向量
    if abs(normal(3)) > 0.5
        tangent1 = [1; 0; 0];
    else
        tangent1 = [0; 0; 1];
    end
    tangent1 = tangent1 - dot(tangent1, normal)*normal;
    tangent1 = tangent1 / norm(tangent1);
    tangent2 = cross(normal, tangent1);
    
    [u,v] = meshgrid(linspace(-size,size,15));
    X = center(1) + u.*tangent1(1) + v.*tangent2(1);
    Y = center(2) + u.*tangent1(2) + v.*tangent2(2);
    Z = center(3) + u.*tangent1(3) + v.*tangent2(3);
    
    surf(X,Y,Z, 'FaceColor',color, 'FaceAlpha',alpha, 'EdgeColor','none');
end