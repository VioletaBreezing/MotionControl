%% 轨迹估计与3D可视化程序
% 功能：根据设定的轨迹，实时可视化真实轨迹、估计轨迹和3D刚体运动
% 基于 cal_eddy_3dof.m 和 simulate_eddy.m 的功能整合，保留完整的3D可视化效果

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
T_total = 3;         % 总仿真时间 (s)
N = round(T_total / Ts);
t = (0:N-1) * Ts;

% 设定轨迹类型
trajectory_type = 'sinusoidal'; % 可选: 'sinusoidal', 'circular'

switch trajectory_type
    case 'sinusoidal'
        % 正弦轨迹
        x_true = 1.5e-3 * sin(2*pi*0.8*t);    % X方向: 1.5mm振幅, 0.8Hz
        y_true = 1.2e-3 * cos(2*pi*0.6*t);    % Y方向: 1.2mm振幅, 0.6Hz  
        rz_true = 3e-3 * sin(2*pi*1.0*t);     % Rz方向: 3mrad振幅, 1.0Hz
        
    case 'circular'
        % 圆形轨迹
        radius = 1.5e-3; % 1.5mm半径
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

% 存储3D可视化数据
vertices_history = cell(N, 1);
P_history = zeros(3, 3, N);

fprintf('开始轨迹估计...\n');
for i = 1:N
    % 当前真实位姿
    phi_true = rz_true(i) * [0; 0; 1];
    R_true = rodrigues(phi_true);
    t_true = [x_true(i); y_true(i); 0];
    
    % 计算传感器测量值（光程）
    h_clean = cal_eddy_dis_3dof(R_true, t_true, param);
    
    % 使用估计算法计算位置
    [x_, y_, rz_] = cal_eddy_3dof(h_clean(1), h_clean(2), h_clean(3), param_cal);
    x_est(i) = x_;
    y_est(i) = y_;
    rz_est(i) = rz_;
    
    % 存储3D可视化数据
    vertices_B = [
        -Lx/2,  Lx/2,  Lx/2,     -Lx/2,     -Lx/2,  Lx/2,  Lx/2,    -Lx/2;
        -Ly/2, -Ly/2,  Ly/2*1.4,  Ly/2*1.4, -Ly/2, -Ly/2,  Ly/2*1.4, Ly/2*1.4;
        -Lz/2, -Lz/2, -Lz/2,     -Lz/2,      Lz/2,  Lz/2,  Lz/2,     Lz/2
    ]';
    vertices_W = (R_true * vertices_B')' + repmat(t_true', 8, 1);
    vertices_history{i} = vertices_W;
    
    P0 = R_true * p0_B + t_true;
    P1 = R_true * p1_B + t_true;
    P2 = R_true * p2_B + t_true;
    P_history(:, :, i) = [P0, P1, P2];
    
    if mod(i, max(1, floor(N/10))) == 0
        fprintf('处理进度: %.1f%%\n', i/N*100);
    end
end
fprintf('轨迹估计完成！\n');

%% ========== 5. 结果可视化 ==========
% 创建主图形窗口
figure('Name', '轨迹估计与3D可视化', 'Position', [100, 100, 1400, 800]);

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

% 子图5: 估计误差
subplot(2, 3, 5);
err_x = (x_true - x_est)*1000;
err_y = (y_true - y_est)*1000;
err_rz = (rz_true - rz_est)*1000;
plot(t, err_x, 'r-', 'LineWidth', 1); hold on;
plot(t, err_y, 'g-', 'LineWidth', 1);
plot(t, err_rz, 'b-', 'LineWidth', 1);
xlabel('时间 (s)'); ylabel('误差 (mm/mrad)');
title('估计误差');
legend('X误差', 'Y误差', 'Rz误差', 'Location', 'best');
grid on;

% 子图6: 3D静态可视化（显示中间时刻的3D场景）
subplot(2, 3, 6);
hold on; axis equal;
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('3D刚体运动可视化 (中间时刻)');

% 选择中间时刻进行3D可视化
mid_idx = floor(N/2);
current_vertices = vertices_history{mid_idx};
current_P = P_history(:, :, mid_idx);

% 定义颜色
colors = lines(3);

% --- 绘制镜面位置 ---
mirror_size = 0.07;
for j = 1:3
    Pi = current_P(:,j);
    % 计算当前镜面法向用于绘制方向指示（简单起见用散点表示中心）
    scatter3(Pi(1), Pi(2), Pi(3), 80, 'filled', 'MarkerFaceColor', colors(j,:), 'MarkerEdgeColor','k');
end

% --- 绘制完整3D刚体 ---
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
    x = current_vertices(idx,1)';
    y = current_vertices(idx,2)';
    z = current_vertices(idx,3)';
    patch(x, y, z, [0.85 0.85 0.85], 'FaceAlpha', 0.6, 'EdgeColor', [0.3 0.3 0.3]);
end

view(30,25); grid on; box on;

%% ========== 6. 性能统计 ==========
fprintf('\n=== 估计性能统计 ===\n');
fprintf('X方向: 最大误差 = %.3f μm, RMS误差 = %.3f μm\n', ...
    norm(err_x, inf)*1000, sqrt(mean(err_x.^2))*1000);
fprintf('Y方向: 最大误差 = %.3f μm, RMS误差 = %.3f μm\n', ...
    norm(err_y, inf)*1000, sqrt(mean(err_y.^2))*1000);
fprintf('Rz方向: 最大误差 = %.3f μrad, RMS误差 = %.3f μrad\n', ...
    norm(err_rz, inf)*1000, sqrt(mean(err_rz.^2))*1000);

%% ========== 7. 3D动画演示 ==========
figure('Name', '3D轨迹动画', 'Position', [200, 200, 800, 600]);

for i = 1:2:N  % 每2帧更新一次，保证流畅性
    clf; hold on; axis equal;
    xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
    title(sprintf('3D刚体运动 - 时间: %.2f s', t(i)));
    
    % 绘制当前帧的3D刚体
    current_vertices = vertices_history{i};
    current_P = P_history(:, :, i);
    
    % 绘制刚体
    for f = 1:size(faces_idx,1)
        idx = faces_idx(f,:);
        x = current_vertices(idx,1)';
        y = current_vertices(idx,2)';
        z = current_vertices(idx,3)';
        patch(x, y, z, [0.85 0.85 0.85], 'FaceAlpha', 0.6, 'EdgeColor', [0.3 0.3 0.3]);
    end
    
    % 绘制镜面
    for j = 1:3
        Pi = current_P(:,j);
        scatter3(Pi(1), Pi(2), Pi(3), 80, 'filled', 'MarkerFaceColor', colors(j,:), 'MarkerEdgeColor','k');
    end
    
    view(30,25); grid on; box on;
    drawnow;
    pause(0.05); % 调整动画速度
end

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
end

function y = my_atan(x)
    y = x - x^3/3;
end