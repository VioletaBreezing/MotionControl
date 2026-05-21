%% 刚体-反射镜光程仿真（完整3D刚体 + 检波面可视化）
% 功能：
%   - 完整6DoF位姿（当前约束：仅Rz + xy平移）
%   - 绘制3D刚体长方体
%   - 绘制每个光路的检波面（过a_i，法向u_i）
%   - 计算并显示光程 h_i

clear; clc; close all;

%% ========== 1. 真实位姿参数（6DoF，显式设零）==========
theta_x = deg2rad(0);           % 绕X轴旋转（设为0）
theta_y = deg2rad(0);           % 绕Y轴旋转（设为0）
theta_z = deg2rad(0.5);  % 绕Z轴旋转（非零）

tx = 0.005;             % X平移
ty = 0;           % Y平移
tz = 0;                % Z平移（设为0）

t = [tx; ty; tz];

%% ========== 2. 构建完整旋转矩阵 ==========
Rx = [1, 0, 0;
      0, cos(theta_x), -sin(theta_x);
      0, sin(theta_x), cos(theta_x)];

Ry = [cos(theta_y), 0, sin(theta_y);
      0, 1, 0;
      -sin(theta_y), 0, cos(theta_y)];

Rz = [cos(theta_z), -sin(theta_z), 0;
      sin(theta_z), cos(theta_z), 0;
      0, 0, 1];

R = Rz * Ry * Rx;  % 固定轴顺序：先Rx，再Ry，再Rz

%% ========== 3. 刚体与镜面定义 ==========
Lx = 0.2; Ly = 0.2; Lz = 0.1;  % 长宽高
d = 0.1;                       % M1-M2间距

% 刚体本体系8个顶点（原点在几何中心）
% 顺序：下底面(z=-Lz/2) -> 上底面(z=+Lz/2)
vertices_B = [
    -Lx/2,  Lx/2,  Lx/2, -Lx/2, -Lx/2,  Lx/2,  Lx/2, -Lx/2;
    -Ly/2, -Ly/2,  Ly/2,  Ly/2, -Ly/2, -Ly/2,  Ly/2,  Ly/2;
    -Lz/2, -Lz/2, -Lz/2, -Lz/2,  Lz/2,  Lz/2,  Lz/2,  Lz/2
]';

% 镜面参考点（位于上表面 z = +Lz/2）
p0_B = [-Lx/2;  0;      0];   % M0: X负向面中心（上边中点）
p1_B = [-d/2;   -Ly/2;  0];   % M1: Y负向面左侧（上边）
p2_B = [ d/2;   -Ly/2;  0];   % M2: Y负向面右侧（上边）

% 镜面法向（指向外）
n0_B = [-1;  0; 0]; n0_B = n0_B / norm(n0_B);
n1_B = [ 0; -1; 0]; n1_B = n1_B / norm(n1_B);
n2_B = [ 0; -1; 0]; n2_B = n2_B / norm(n2_B);

% 转换到世界系
vertices_W = (R * vertices_B')' + repmat(t', 8, 1);
P0 = R * p0_B + t;  N0 = R * n0_B;
P1 = R * p1_B + t;  N1 = R * n1_B;
P2 = R * p2_B + t;  N2 = R * n2_B;

%% ========== 4. 入射光线与检波面 ==========
a0 = [-0.3;  0.0;  0.0];
u0 = [ 1.0;  0.0;  0.0]; u0 = u0 / norm(u0);

a1 = [-d/2; -0.3; 0.0];
a2 = [ d/2; -0.3; 0.0];
u1 = [0.0;  1.0;  0.0]; u1 = u1 / norm(u1);
u2 = [0.0;  1.0;  0.0]; u2 = u2 / norm(u2);

rays.a = {a0, a1, a2};
rays.u = {u0, u1, u2};

%% ========== 5. 光程计算 ==========
h = zeros(3,1);
h_ = zeros(3,1);
c_points = cell(3,1);
b_points = cell(3,1);

h_(1) = cal_light_dis(R, t, a0, u0, p0_B, n0_B);
h_(2) = cal_light_dis(R, t, a1, u1, p1_B, n1_B);
h_(3) = cal_light_dis(R, t, a2, u2, p2_B, n2_B);

for i = 1:3
    Pi = eval(['P' num2str(i-1)]);  % P0, P1, P2
    Ni = eval(['N' num2str(i-1)]);
    ai = rays.a{i};
    ui = rays.u{i};
    
    Ni = Ni / norm(Ni);
    beta = dot(ui, Ni);
    if abs(beta) < 1e-8
        error('光路 %d: 入射光与镜面平行！', i-1);
    end
    
    si = dot(Pi - ai, Ni) / beta;
    ci = ai + si * ui;
    
    ui_prime = ui - 2 * beta * Ni;
    ui_prime = ui_prime / norm(ui_prime);
    
    denom = dot(ui_prime, ui);
    if abs(denom) < 1e-8
        warning('光路 %d: 反射光与检波面平行', i-1);
        ui_param = 0;
    else
        ui_param = -si / denom;
    end
    bi = ci + ui_param * ui_prime;
    
    h(i) = (si + ui_param) / 2;
    
    c_points{i} = ci;
    b_points{i} = bi;
end

%% ========== 6. 可视化 ==========
figure('Color','w'); hold on; axis equal;
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('刚体-反射镜系统仿真（含检波面）');

% --- 绘制镜面（小方块，贴在刚体表面）---
mirror_size = 0.07;
mirror_centers = [P0, P1, P2];
mirror_normals = [N0, N1, N2];
for i = 1:3
    Pi = mirror_centers(:,i);
    Ni = mirror_normals(:,i);
    if abs(Ni(3)) < 0.9
        v1 = cross(Ni, [0;0;1]);
    else
        v1 = cross(Ni, [1;0;0]);
    end
    v1 = v1 / norm(v1) * mirror_size/2;
    v2 = cross(Ni, v1);
    v2 = v2 / norm(v2) * mirror_size/2;
    mc = Pi + [v1+v2, -v1+v2, -v1-v2, v1-v2];
    patch(mc(1,:), mc(2,:), mc(3,:), [1 0.5 0.5], ...
          'EdgeColor','r', 'LineWidth',1.5, 'FaceAlpha',0.7);
end

% --- 绘制检波面（每个光路一个）---
colors = lines(3);
plane_size = 0.07;  % 检波面尺寸
for i = 1:3
    ai = rays.a{i};
    ui = rays.u{i};
    
    % 构造检波面局部坐标系（两个正交于ui的向量）
    if abs(ui(3)) < 0.9
        e1 = cross(ui, [0;0;1]);
    else
        e1 = cross(ui, [1;0;0]);
    end
    e1 = e1 / norm(e1) * plane_size/2;
    e2 = cross(ui, e1);
    e2 = e2 / norm(e2) * plane_size/2;
    
    % 检波面四个角点
    plane_corners = ai + [e1+e2, -e1+e2, -e1-e2, e1-e2];
    
    % 绘制半透明平面
    patch(plane_corners(1,:), plane_corners(2,:), plane_corners(3,:), ...
          colors(i,:), 'FaceAlpha', 0.2, 'EdgeColor', colors(i,:), 'LineStyle', ':');
    
    % 标注"检波面 i"
    text(ai(1)+0.02, ai(2)+0.02, ai(3)+0.02, sprintf('Detector %d', i-1), ...
         'Color', colors(i,:), 'FontSize', 9, 'FontWeight','bold');
end

% --- 绘制光线路径 ---
for i = 1:3
    ai = rays.a{i};
    ci = c_points{i};
    bi = b_points{i};
    col = colors(i,:);
    
    % 入射光（虚线）
    plot3([ai(1), ci(1)], [ai(2), ci(2)], [ai(3), ci(3)], '--', 'Color', col, 'LineWidth',1.5);
    % 反射光（实线）
    plot3([ci(1), bi(1)], [ci(2), bi(2)], [ci(3), bi(3)], '-', 'Color', col, 'LineWidth',1.5);
    
    % 关键点
    scatter3(ai(1), ai(2), ai(3), 60, 'filled', 'MarkerFaceColor', col, 'MarkerEdgeColor','k');
    scatter3(ci(1), ci(2), ci(3), 60, 'ko', 'MarkerFaceColor', 'g', 'MarkerEdgeColor','k');
    scatter3(bi(1), bi(2), bi(3), 60, 's', 'MarkerFaceColor', col, 'MarkerEdgeColor','k');
end

% --- 绘制完整3D刚体 ---
% 定义6个面的顶点索引（每行一个面，4个顶点）
faces_idx = [
    1 2 3 4;  % 下底面
    5 6 7 8;  % 上底面
    1 2 6 5;  % 前面 (Y负)
    3 4 8 7;  % 后面 (Y正)
    1 4 8 5;  % 左面 (X负)
    2 3 7 6   % 右面 (X正)
];

for f = 1:size(faces_idx,1)
    idx = faces_idx(f,:);
    x = vertices_W(idx,1)';
    y = vertices_W(idx,2)';
    z = vertices_W(idx,3)';
    patch(x, y, z, [0.85 0.85 0.85], 'FaceAlpha', 0.4, 'EdgeColor', [0.3 0.3 0.3]);
end

legend_str = {'M_0','M_1','M_2'};
for i=1:3
    legend_str{end+1} = sprintf('检波面 %d', i-1);
end

legend_str{end+1} = "入射光";
legend_str{end+1} = "反射光";

legend(legend_str, 'Location','bestoutside');

view(30,25); grid on; box on;

%% ========== 7. 输出结果 ==========
fprintf('\n=== 位姿参数 ===\n');
fprintf('旋转: Rx=%.2f°, Ry=%.2f°, Rz=%.2f°\n', rad2deg([theta_x, theta_y, theta_z]));
fprintf('平移: tx=%.1f mm, ty=%.1f mm, tz=%.1f mm\n', tx*1000, ty*1000, tz*1000);
fprintf('\n光程测量值:\n');
for i=1:3
    fprintf('  h_%d = %.10f m\n', i-1, h(i));
end

function h = cal_light_dis(R, t, a, u, p, n)
    N = R * n;
    beta = dot(u, N);
    ga = dot(p, n);
    d = t - a;
    de = dot(d, N);
    f = beta / (2*beta^2 - 1);
    g = ga + de;
    h = f * g;
end