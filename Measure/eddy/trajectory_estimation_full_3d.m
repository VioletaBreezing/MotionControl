%% 轨迹估计与完整3D可视化程序
% 功能：根据设定的轨迹，一边可视化（包含镜面、检波面、光线路径），一边估计轨迹
% 完全保留 simulate_eddy.m 中的3D可视化效果

clear; clc; close all;

%% ========== 1. 系统参数设置 ==========
% 电涡流测量系统参数
param_cal.ax = 0.09;
param_cal.ay1 = 0.107; 
param_cal.ay2 = 0.107;
param_cal.bx = 0.1975;
param_cal.by = 0.116;

% 刚体几何参数
ax = param_cal.ax; ay1 = param_cal.ay1; ay2 = param_cal.ay2;
bx = param_cal.bx; by = param_cal.by;
Lx = 2*bx; Ly = 2*by; Lz = 0.1;  % 长宽高
d = ay1 + ay2;                   % M1-M2间距

% 轨迹参数
Ts = 1e-2;           % 采样时间 (s)
T_total = 8;         % 总仿真时间 (s)
N = round(T_total / Ts);
t = (0:N-1) * Ts;

% 设定轨迹类型
trajectory_type = 'sinusoidal'; % 可选: 'sinusoidal', 'circular'
noise_amp = 3e-10/6;

switch trajectory_type
    case 'sinusoidal'
        % 正弦轨迹
        x_true = 1.5e-3 * sin(2*pi*0.8*t);    % X方向: 1.5mm振幅, 0.8Hz
        y_true = 1.2e-3 * sin(2*pi*0.6*t);    % Y方向: 1.2mm振幅, 0.6Hz  
        rz_true = 3e-3 * sin(2*pi*1.0*t);     % Rz方向: 3mrad振幅, 1.0Hz
        
    case 'circular'
        % 圆形轨迹
        radius = 4e-3; % 1.5mm半径
        omega = 2*pi*0.5; % 0.5Hz
        x_true = radius * cos(omega * t);
        y_true = radius * sin(omega * t);
        rz_true = 2e-3 * sin(omega * t);
end

%% ========== 2. 传感器几何配置 ==========
% 镜面参考点（本体系）
p0_B = [Lx/2;  ax;      0];   % M0: X负向面中心（上边中点）
p1_B = [-d/2;   -Ly/2;  0];   % M1: Y负向面左侧（上边）
p2_B = [ d/2;   -Ly/2;  0];   % M2: Y负向面右侧（上边）

% 镜面法向（指向外）
n0_B = [ 1;  0; 0]; n0_B = n0_B / norm(n0_B);
n1_B = [ 0; -1; 0]; n1_B = n1_B / norm(n1_B);
n2_B = [ 0; -1; 0]; n2_B = n2_B / norm(n2_B);

% 入射光线与检波面
a0 = [ 1.4*bx;  ax;  0.0];
u0 = [-1.0;  0;  0]; u0 = u0 / norm(u0);

a1 = [-ay1; -by*1.4; 0.0];
a2 = [ ay2; -by*1.4; 0.0];
u1 = [0;  1.0;  0]; u1 = u1 / norm(u1);
u2 = [0;  1.0;  0]; u2 = u2 / norm(u2);

rays.a = {a0, a1, a2};
rays.u = {u0, u1, u2};

% 构建参数结构体用于估计
param.a = {a0, a1, a2};
param.u = {u0, u1, u2};
param.p = {p0_B, p1_B, p2_B};
param.n = {n0_B, n1_B, n2_B};

%% ========== 3. 计算死程（零位时光程）==========
R_eye = eye(3,3);
t_zero = zeros(3,1);
l0 = zeros(3,1);
l0(1) = cal_eddy_len(R_eye, t_zero, a0, u0, p0_B, n0_B);
l0(2) = cal_eddy_len(R_eye, t_zero, a1, u1, p1_B, n1_B);
l0(3) = cal_eddy_len(R_eye, t_zero, a2, u2, p2_B, n2_B);
param.l0 = l0;

%% ========== 4. 轨迹估计 ==========
x_est = zeros(1, N);
y_est = zeros(1, N);
rz_est = zeros(1, N);

% 存储3D可视化数据（选择中间时刻进行可视化）
mid_idx = floor(N/2);
phi_mid = rz_true(mid_idx) * [0; 0; 1];
R_mid = rodrigues(phi_mid);
t_mid = [x_true(mid_idx); y_true(mid_idx); 0];

fprintf('开始轨迹估计...\n');
for i = 1:N
    % 当前真实位姿
    phi_true = rz_true(i) * [0; 0; 1];
    R_true = rodrigues(phi_true);
    t_true = [x_true(i); y_true(i); 0];
    
    % 计算传感器测量值（光程）
    h_clean = cal_eddy_dis_3dof(R_true, t_true, param);

    % 测量噪声
    h_meas = h_clean +  + noise_amp * randn(3, 1);
    
    % 使用估计算法计算位置
    [x_, y_, rz_] = cal_eddy_3dof(h_meas(1), h_meas(2), h_meas(3), struct('ax', ax, 'ay1', ay1, 'ay2', ay2, 'bx', bx, 'by', by));
    x_est(i) = x_;
    y_est(i) = y_;
    rz_est(i) = rz_;
    
    if mod(i, max(1, floor(N/10))) == 0
        fprintf('处理进度: %.1f%%\n', i/N*100);
    end
end
fprintf('轨迹估计完成！\n');

%% ========== 5. 结果可视化 ==========
% 创建主图形窗口
figure('Name', '轨迹估计与完整3D可视化', 'Position', [100, 100, 1400, 800]);

% 子图1: X方向轨迹对比
subplot(2, 3, 1);
plot(t, x_true*1000, 'b-', 'LineWidth', 1.5); hold on;
plot(t, x_est*1000, 'r--', 'LineWidth', 1.5);
xlabel('时间 (s)'); ylabel('X (mm)');
title('X方向轨迹');
legend('真实', '估计', 'Location', 'best');
grid on;

% 子图2: Y方向轨迹对比
subplot(2, 3, 2);
plot(t, y_true*1000, 'b-', 'LineWidth', 1.5); hold on;
plot(t, y_est*1000, 'r--', 'LineWidth', 1.5);
xlabel('时间 (s)'); ylabel('Y (mm)');
title('Y方向轨迹');
legend('真实', '估计', 'Location', 'best');
grid on;

% 子图3: Rz方向轨迹对比
subplot(2, 3, 3);
plot(t, rz_true*1000, 'b-', 'LineWidth', 1.5); hold on;
plot(t, rz_est*1000, 'r--', 'LineWidth', 1.5);
xlabel('时间 (s)'); ylabel('Rz (mrad)');
title('Rz方向轨迹');
legend('真实', '估计', 'Location', 'best');
grid on;

% 子图4: XY平面轨迹
subplot(2, 3, 4);
plot(x_true*1000, y_true*1000, 'b-', 'LineWidth', 1.5); hold on;
plot(x_est*1000, y_est*1000, 'r--', 'LineWidth', 1.5);
xlabel('X (mm)'); ylabel('Y (mm)');
title('XY平面轨迹');
legend('真实', '估计', 'Location', 'best');
axis equal; grid on;

% 子图5: 估计误差（单位：nm）
subplot(2, 3, 5);
err_x_nm = (x_true - x_est) * 1e9;  % 转换为nm
err_y_nm = (y_true - y_est) * 1e9;  % 转换为nm
err_rz_nrad = (rz_true - rz_est) * 1e9;  % 转换为nrad
plot(t, err_x_nm, 'r-', 'LineWidth', 1); hold on;
plot(t, err_y_nm, 'g-', 'LineWidth', 1);
plot(t, err_rz_nrad, 'b-', 'LineWidth', 1);
xlabel('时间 (s)'); ylabel('误差 (nm/nrad)');
title('估计误差');
legend('X误差 (nm)', 'Y误差 (nm)', 'Rz误差 (nrad)', 'Location', 'best');
grid on;

%% ========== 6. 完整3D可视化（像simulate_eddy.m一样）==========
% 选择最终时刻进行完整3D可视化
final_idx = N;
phi_final = rz_true(final_idx) * [0; 0; 1];
R_final = rodrigues(phi_final);
t_final = [x_true(final_idx); y_true(final_idx); 0];

% 计算最终时刻的3D数据
vertices_B = [
    -Lx/2,  Lx/2,  Lx/2,     -Lx/2,     -Lx/2,  Lx/2,  Lx/2,    -Lx/2;
    -Ly/2, -Ly/2,  Ly/2*1.4,  Ly/2*1.4, -Ly/2, -Ly/2,  Ly/2*1.4, Ly/2*1.4;
    -Lz/2, -Lz/2, -Lz/2,     -Lz/2,      Lz/2,  Lz/2,  Lz/2,     Lz/2
]';
vertices_W = (R_final * vertices_B')' + repmat(t_final', 8, 1);

P0 = R_final * p0_B + t_final;  N0 = R_final * n0_B;
P1 = R_final * p1_B + t_final;  N1 = R_final * n1_B;
P2 = R_final * p2_B + t_final;  N2 = R_final * n2_B;

% 子图6: 完整3D可视化
subplot(2, 3, 6);
hold on; axis equal;
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('完整3D可视化（最终时刻）');

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
c_points = cell(3,1);
b_points = cell(3,1);
l = zeros(3,1);
h = zeros(3,1);

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
    
    l(i) = si;  % 电涡流的
    h(i) = l(i) - l0(i);
    if i == 1
        h(i) = -h(i);
    end
    
    c_points{i} = ci;
    b_points{i} = bi;
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

%% ========== 7. 性能统计 ==========
fprintf('\n=== 估计性能统计 ===\n');
fprintf('X方向: 最大误差 = %.3e nm, RMS误差 = %.3e nm\n', ...
    norm(err_x_nm, inf), sqrt(mean(err_x_nm.^2)));
fprintf('Y方向: 最大误差 = %.3e nm, RMS误差 = %.3e nm\n', ...
    norm(err_y_nm, inf), sqrt(mean(err_y_nm.^2)));
fprintf('Rz方向: 最大误差 = %.3e nrad, RMS误差 = %.3e nrad\n', ...
    norm(err_rz_nrad, inf), sqrt(mean(err_rz_nrad.^2)));

%% ========== 辅助函数定义 ==========
function h = cal_eddy_dis_3dof(R, t, param)
    h = zeros(3, 1);
    for i = 1:3
        a = param.a{i};
        u = param.u{i};
        p = param.p{i};
        n = param.n{i};
        l0 = param.l0(i);
        h(i) = cal_eddy_dis(R, t, a, u, p, n, l0);
    end
    h(1) = -h(1);
end

function l = cal_eddy_len(R, t, a, u, p, n)
    N = R * n;
    P = R * p + t;
    beta = dot(u, N);
    l = dot((P - a), N) / beta;
end

function h = cal_eddy_dis(R, t, a, u, p, n, l0)
    l = cal_eddy_len(R, t, a, u, p, n);
    h = l - l0;
end

function R = rodrigues(phi)
    theta = norm(phi);
    if theta < 1e-8
        R = eye(3);
    else
        k = phi / theta;
        K = skew_matrix(k);
        R = cos(theta)*eye(3) + (1-cos(theta))*(k*k') + sin(theta)*K;
    end
end

function maa = skew_matrix(a)
    maa = [  0, -a(3),  a(2);
          a(3),   0, -a(1);
         -a(2), a(1),   0];
end

function [x, y, rz] = cal_eddy_3dof(hx, hy1, hy2, param)
    ax = param.ax; ay1 = param.ay1; ay2 = param.ay2;
    bx = param.bx; by = param.by;
    ay = ay1 + ay2;

    x = hx;
    y = ay2/ay * hy1 + ay1/ay * hy2;
    rz = my_atan((hy2 - hy1) / ay);

    for iter = 1:20
        armx  = -y + ax  - bx * my_tan(rz/2);
        army1 =  x + ay1 + by * my_tan(rz/2);
        army2 = -x + ay2 - by * my_tan(rz/2);
    
        abbey_errx  = armx  * my_tan(rz);
        abbey_erry1 = army1 * my_tan(rz);
        abbey_erry2 = army2 * my_tan(rz);
    
        x_ = hx + abbey_errx;
        y_ = ay2/ay * (hy1 + abbey_erry1) + ay1/ay * (hy2 - abbey_erry2);

        v  = [x, y];
        dv = [x_-x, y_-y];

        x = x_;
        y = y_;

        if norm(dv) / norm(v) < 1e-6
            break;
        end
    end
    if iter == 20
        fprintf("iter = %d, cannot converge!!", iter);
    end
end

function y = my_tan(x)
    y = x + x^3/3;
    y = tan(x);
end

function y = my_atan(x)
    y = x - x^3/3;
    y = atan(x);
end