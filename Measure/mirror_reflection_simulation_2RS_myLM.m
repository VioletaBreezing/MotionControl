%% 刚体-反射镜光程仿真（完整3D刚体 + 检波面可视化）
% 功能：
%   - 完整6DoF位姿（当前约束：仅Rz + xy平移）
%   - 绘制3D刚体长方体
%   - 绘制每个光路的检波面（过a_i，法向u_i）
%   - 计算并显示光程 h_i

clear; clc; close all;

path_include;

%% ========== 1. 真实位姿参数（6DoF，显式设零）==========
theta_x = deg2rad(0);           % 绕X轴旋转（设为0）
theta_y = deg2rad(0);           % 绕Y轴旋转（设为0）
theta_z = deg2rad(20e-6/pi*180);         % 绕Z轴旋转（非零）

tx = 1e-3;             % X平移
ty = -1e-3;           % Y平移
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
u0 = [ 1.0;  0.01;  -0.01]; u0 = u0 / norm(u0);

a1 = [-d/2; -0.3; 0.0];
a2 = [ d/2; -0.3; 0.0];
u1 = [-0.01;  1.0;  0.02]; u1 = u1 / norm(u1);
u2 = [0.04;  1.0;  0.03]; u2 = u2 / norm(u2);

rays.a = {a0, a1, a2};
rays.u = {u0, u1, u2};

%% ========== 5. 光程计算 ==========
h = zeros(3,1);
h_ = zeros(3,1);
jac= zeros(3,6);
c_points = cell(3,1);
b_points = cell(3,1);

h_(1) = cal_light_dis(R, t, a0, u0, p0_B, n0_B);
h_(2) = cal_light_dis(R, t, a1, u1, p1_B, n1_B);
h_(3) = cal_light_dis(R, t, a2, u2, p2_B, n2_B);

jac(1,:) = calc_jacobian(R, t, a0, u0, p0_B, n0_B);
jac(2,:) = calc_jacobian(R, t, a1, u1, p1_B, n1_B);
jac(3,:) = calc_jacobian(R, t, a2, u2, p2_B, n2_B);

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
% close all;
%% ========== 8. LM 优化：从 h_ 反演位姿 ==========
fprintf('\n=== 开始 LM 优化 ===\n');

% --- 已知配置（来自前面定义）---
p_B = {p0_B, p1_B, p2_B};
n_B = {n0_B, n1_B, n2_B};
a_list = [a0, a1, a2];   % 3×3
u_list = [u0, u1, u2];   % 3×3
h_meas = h_;             % 3×1，真实测量值

% --- 初始猜测（故意设错，测试收敛性）---
t0 = [0; 0; 0];                     % 初始平移 = 0
phi0 = [0; 0; deg2rad(0)];         % 初始旋转 = 0（比真实值小）
xi = [t0; phi0];                   % 6×1: [tx; ty; tz; phix; phiy; phiz]

% --- LM 参数 ---
max_iter = 5;
lambda = 1e-3;
nu = 2;
TolX = 1e-9;
TolR = 1e-9;
prev_cost = inf;

% --- 预计算常量（镜面在本体系中的 gamma = p' * n）---
gamma_list = zeros(3,1);
for i = 1:3
    gamma_list(i) = p_B{i}' * n_B{i};
end

% 当前真实位姿
t_true = [tx; ty; tz];
phi_true = R_to_phi(R);
R_true = R;

args.a = {a0, a1, a2};
args.u = {u0, u1, u2};
args.p = {p0_B, p1_B, p2_B};
args.n = {n0_B, n1_B, n2_B};
args.h = h_;
args.OutDim = 3;

N = 200;
tic;  % ←←← 开始计时

for i = 1:N
[x_hat, resnorm, converged, iter] = Optimizer_LM(@Get_RS2_IFM_MS_ResFcn, args, ...
             @Get_RS2_IFM_MS_Jacobian, args, ...
             zeros(6,1), TolX, TolR, 1e10,...
             max_iter, false);
end

elapsed_time_us = toc * 1000000 / N;  % ←←← 结束计时，转为微秒平均耗时
fprintf('%d 次优化平均耗时: %.3f us\n', N, elapsed_time_us);

% --- 提取最终结果 ---
T_opt = se3_to_SE3(x_hat);
% t_opt = x_hat(1:3);
phi_opt = x_hat(4:6);
% R_opt = rodrigues(phi_opt);
t_opt = T_opt(1:3, 4);
R_opt = T_opt(1:3, 1:3);

true_t = t;
true_phi = R_to_phi(R_true);
trans_error_m = (t_opt - t);                     % 平移误差 → m
rot_error_rad  = phi_opt - true_phi;             % 旋转李代数误差 → 弧度

% 转换为欧拉角（Z-Y-X）
[theta_x_opt, theta_y_opt, theta_z_opt] = rot2euler(R_opt);

fprintf('\n=== 优化结果 vs 真实值 ===\n');
fprintf('平移 (mm):     真实 [%.3f, %.3f, %.3f] -> 估计 [%.3f, %.3f, %.3f]\n', ...
        tx*1000, ty*1000, tz*1000, t_opt(1)*1000, t_opt(2)*1000, t_opt(3)*1000);
fprintf('旋转 (deg):    真实 [%.3f, %.3f, %.3f] -> 估计 [%.3f, %.3f, %.3f]\n', ...
        rad2deg([theta_x, theta_y, theta_z]), ...
        rad2deg([theta_x_opt, theta_y_opt, theta_z_opt]));
fprintf('平移误差 (nm): %.2e, %.2e, %.2e\n', trans_error_m * 1e9);
fprintf('so(3)误差 (nrad): %.2e, %.2e, %.2e\n', rot_error_rad * 1e9);

%% 函数定义

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

function jac = calc_jacobian(R, t, a, u, p, n) 
    N = R * n;
    beta = dot(u, N);
    ga = dot(p, n);
    d = t - a;
    de = dot(d, N);
    f = beta / (2*beta^2 - 1);
    g = ga + de;
    alph = -(2*beta^2 + 1) / (2*beta^2 - 1)^2;

    pt = f * N';
    ph = -transpose(alph * g * u + f * d) * skew_matrix(N);

    jac = [pt, ph];
end

function maa = skew_matrix(a)
    maa = [  0, -a(3),  a(2);
          a(3),   0, -a(1);
         -a(2), a(1),   0];
end

% --- 辅助函数：从旋转矩阵提取 ZYX 欧拉角 ---
function [rx, ry, rz] = rot2euler(R)
    % R = Rz * Ry * Rx
    ry = atan2(-R(3,1), sqrt(R(1,1)^2 + R(2,1)^2));
    if abs(ry - pi/2) < 1e-6
        rx = 0;
        rz = atan2(R(1,2), R(1,3));
    elseif abs(ry + pi/2) < 1e-6
        rx = 0;
        rz = atan2(R(1,2), R(1,3));
    else
        rx = atan2(R(3,2), R(3,3));
        rz = atan2(R(2,1), R(1,1));
    end
end

function R = rodrigues(phi)
    theta = norm(phi);
    if theta < 1e-8
        R = eye(3);
    else
        k = phi / theta;
        K = skew_matrix(k);
        % R = eye(3) + sin(theta)*K + (1-cos(theta))*K*K;
        R = cos(theta)*eye(3) + (1-cos(theta))*(k*k') + sin(theta)*K;
    end
end

function phi = R_to_phi(R)
% 将 SO(3) 中的旋转矩阵 R 转换为 so(3) 中的旋转向量 phi
% 输入: R (3x3 旋转矩阵)
% 输出: phi (3x1 向量，单位：弧度)

    % 确保 R 是有效的旋转矩阵（可选）
    assert(all(size(R) == [3,3]), 'R must be 3x3');
    
    % 计算旋转角 theta
    cos_theta = (trace(R) - 1) / 2;
    cos_theta = max(-1.0, min(1.0, cos_theta)); % 防止数值误差导致 acos 越界
    theta = acos(cos_theta);

    if theta < 1e-12
        % 情况1: 无旋转 → phi = [0;0;0]
        phi = zeros(3,1);
    elseif theta > pi - 1e-6
        % 情况2: 接近180度旋转（数值不稳定，需特殊处理）
        % 利用 R + I 的最大对角元确定主轴
        A = (R + eye(3)) / 2;
        [~, idx] = max(diag(A));
        n = A(:, idx) / norm(A(:, idx));
        phi = pi * n;
    else
        % 情况3: 一般情况
        sin_theta = sin(theta);
        nx = (R(3,2) - R(2,3)) / (2*sin_theta);
        ny = (R(1,3) - R(3,1)) / (2*sin_theta);
        nz = (R(2,1) - R(1,2)) / (2*sin_theta);
        n = [nx; ny; nz];
        phi = theta * n;
    end
end

function sovle_linear_equation_cholesky()
end

function J = Get_RS2_IFM_MS_Jacobian(xi, args)
    T = se3_to_SE3(xi);

    R = T(1:3, 1:3);
    t = T(1:3, 4);
    J = zeros(args.OutDim, 6);
    
    for i = 1:args.OutDim
        u = args.u{i};
        a = args.a{i};
        n = args.n{i};
        p = args.p{i};
        J(i, :) = calc_jacobian(R, t, a, u, p, n);
        J(i, 3) = 0;
        J(i, 4) = 0;
        J(i, 5) = 0;
    end
end

function residul = Get_RS2_IFM_MS_ResFcn(xi, args)
    T = se3_to_SE3(xi);
    R = T(1:3, 1:3);
    t = T(1:3, 4);
    % R = args.R;
    % t = args.t;
    h_ = zeros(args.OutDim, 1);

    for i = 1:args.OutDim
        u = args.u{i};
        a = args.a{i};
        n = args.n{i};
        p = args.p{i};
        h_(i) = cal_light_dis(R, t, a, u, p, n);
    end

    residul = h_ - args.h;
end