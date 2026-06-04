%% 轨迹估计与可视化程序
% 功能：根据设定的轨迹，实时可视化真实轨迹和估计轨迹
% 基于 cal_eddy_3dof.m 和 simulate_eddy.m 的功能整合

clear; clc; close all;

%% ========== 1. 系统参数设置 ==========
% 电涡流测量系统参数
param_cal.ax = 0.09;
param_cal.ay1 = 0.107; 
param_cal.ay2 = 0.107;
param_cal.bx = 0.1975;
param_cal.by = 0.116;

% 轨迹参数
Ts = 1e-2;           % 采样时间 (s)
T_total = 5;         % 总仿真时间 (s)
N = round(T_total / Ts);
t = (0:N-1) * Ts;

% 设定轨迹类型（可选择不同的轨迹）
trajectory_type = 'sinusoidal'; % 可选: 'sinusoidal', 'circular', 'linear'

switch trajectory_type
    case 'sinusoidal'
        % 正弦轨迹
        x_true = 2e-3 * sin(2*pi*0.5*t);      % X方向: 2mm振幅, 0.5Hz
        y_true = 1.5e-3 * cos(2*pi*0.3*t);    % Y方向: 1.5mm振幅, 0.3Hz  
        rz_true = 5e-3 * sin(2*pi*0.8*t);     % Rz方向: 5mrad振幅, 0.8Hz
        
    case 'circular'
        % 圆形轨迹
        radius = 2e-3; % 2mm半径
        omega = 2*pi*0.5; % 0.5Hz
        x_true = radius * cos(omega * t);
        y_true = radius * sin(omega * t);
        rz_true = 2e-3 * sin(omega * t);
        
    case 'linear'
        % 线性轨迹（带斜坡）
        x_true = 1e-3 * t; % 1mm/s
        y_true = 0.5e-3 * t; % 0.5mm/s
        rz_true = 1e-3 * t; % 1mrad/s
end

%% ========== 2. 初始化传感器几何参数 ==========
% 镜面参考点（本体系）
p0_B = [param_cal.bx; param_cal.ax; 0];   
p1_B = [-(param_cal.ay1 + param_cal.ay2)/2; -param_cal.by; 0];   
p2_B = [(param_cal.ay1 + param_cal.ay2)/2; -param_cal.by; 0];   

% 镜面法向（指向外）
n0_B = [1; 0; 0]; n0_B = n0_B / norm(n0_B);
n1_B = [0; -1; 0]; n1_B = n1_B / norm(n1_B);
n2_B = [0; -1; 0]; n2_B = n2_B / norm(n2_B);

% 入射光线参数
a0 = [1.4*param_cal.bx; param_cal.ax; 0.0];
u0 = [-1.0; 0; 0]; u0 = u0 / norm(u0);

a1 = [-param_cal.ay1; -param_cal.by*1.4; 0.0];
a2 = [param_cal.ay2; -param_cal.by*1.4; 0.0];
u1 = [0; 1.0; 0]; u1 = u1 / norm(u1);
u2 = [0; 1.0; 0]; u2 = u2 / norm(u2);

% 存储参数结构
param.a = {a0, a1, a2};
param.u = {u0, u1, u2};
param.p = {p0_B, p1_B, p2_B};
param.n = {n0_B, n1_B, n2_B};

%% ========== 3. 计算死程（零位光程）==========
R_eye = eye(3,3);
t_zero = zeros(3,1);
l0(1) = cal_eddy_len(R_eye, t_zero, a0, u0, p0_B, n0_B);
l0(2) = cal_eddy_len(R_eye, t_zero, a1, u1, p1_B, n1_B);
l0(3) = cal_eddy_len(R_eye, t_zero, a2, u2, p2_B, n2_B);
param.l0 = l0;

%% ========== 4. 轨迹估计 ==========
x_est = zeros(1, N);  % 初始化为行向量
y_est = zeros(1, N);  % 初始化为行向量
rz_est = zeros(1, N); % 初始化为行向量

fprintf('开始轨迹估计...\n');
for i = 1:N
    % 构建当前时刻的真实位姿
    phi = rz_true(i) * [0; 0; 1];
    R_true = rodrigues(phi);
    t_true = [x_true(i); y_true(i); 0];
    
    % 计算传感器测量值
    h_clean = cal_eddy_dis_3dof(R_true, t_true, param);
    
    % 使用估计算法
    [x_, y_, rz_] = cal_eddy_3dof(h_clean(1), h_clean(2), h_clean(3), param_cal);
    
    x_est(i) = x_;
    y_est(i) = y_;
    rz_est(i) = rz_;
    
    if mod(i, 100) == 0
        fprintf('处理进度: %.1f%%\n', i/N*100);
    end
end
fprintf('轨迹估计完成！\n');

%% ========== 5. 结果可视化 ==========
% 所有向量现在都是行向量，无需额外转换

% 创建主图形窗口
figure('Name', '轨迹估计与可视化', 'Position', [100, 100, 1200, 800]);

% 子图1: X方向轨迹
subplot(3, 2, 1);
plot(t, x_true*1000, 'b-', 'LineWidth', 2); hold on;
plot(t, x_est*1000, 'r--', 'LineWidth', 2);
xlabel('时间 (s)'); ylabel('X位置 (mm)');
title('X方向轨迹');
legend('真实值', '估计值', 'Location', 'best');
grid on;

% 子图2: Y方向轨迹  
subplot(3, 2, 3);
plot(t, y_true*1000, 'b-', 'LineWidth', 2); hold on;
plot(t, y_est*1000, 'r--', 'LineWidth', 2);
xlabel('时间 (s)'); ylabel('Y位置 (mm)');
title('Y方向轨迹');
legend('真实值', '估计值', 'Location', 'best');
grid on;

% 子图3: Rz方向轨迹
subplot(3, 2, 5);
plot(t, rz_true*1000, 'b-', 'LineWidth', 2); hold on;
plot(t, rz_est*1000, 'r--', 'LineWidth', 2);
xlabel('时间 (s)'); ylabel('Rz角度 (mrad)');
title('Rz方向轨迹');
legend('真实值', '估计值', 'Location', 'best');
grid on;

% 子图4: XY平面轨迹
subplot(3, 2, 2);
plot(x_true*1000, y_true*1000, 'b-', 'LineWidth', 2); hold on;
plot(x_est*1000, y_est*1000, 'r--', 'LineWidth', 2);
xlabel('X位置 (mm)'); ylabel('Y位置 (mm)');
title('XY平面轨迹');
legend('真实值', '估计值', 'Location', 'best');
axis equal; grid on;

% 子图5: 估计误差
subplot(3, 2, 4);
plot(t, (x_true - x_est)*1000, 'r-', 'LineWidth', 1.5); hold on;
plot(t, (y_true - y_est)*1000, 'g-', 'LineWidth', 1.5);
plot(t, (rz_true - rz_est)*1000, 'b-', 'LineWidth', 1.5);
xlabel('时间 (s)'); ylabel('误差 (mm/mrad)');
title('估计误差');
legend('X误差', 'Y误差', 'Rz误差', 'Location', 'best');
grid on;

% 子图6: 3D轨迹动画（显示最后几个点）
subplot(3, 2, 6);
plot3(x_true(end-50:end)*1000, y_true(end-50:end)*1000, rz_true(end-50:end)*1000, 'b-', 'LineWidth', 2); hold on;
plot3(x_est(end-50:end)*1000, y_est(end-50:end)*1000, rz_est(end-50:end)*1000, 'r--', 'LineWidth', 2);
xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Rz (mrad)');
title('3D轨迹（最近50个点）');
legend('真实值', '估计值', 'Location', 'best');
grid on;

%% ========== 6. 性能统计 ==========
err_x = x_true - x_est;
err_y = y_true - y_est;
err_rz = rz_true - rz_est;

fprintf('\n=== 估计性能统计 ===\n');
fprintf('X方向: RMSE = %.3e nm, 最大误差 = %.3e nm\n', ...
    sqrt(mean(err_x.^2))*1e9, max(abs(err_x))*1e9);
fprintf('Y方向: RMSE = %.3e nm, 最大误差 = %.3e nm\n', ...
    sqrt(mean(err_y.^2))*1e9, max(abs(err_y))*1e9);
fprintf('Rz方向: RMSE = %.3e nrad, 最大误差 = %.3e nrad\n', ...
    sqrt(mean(err_rz.^2))*1e9, max(abs(err_rz))*1e9);

%% ========== 7. 实时动画演示（可选）==========
show_animation = true;
if show_animation
    figure('Name', '实时轨迹动画', 'Position', [200, 200, 800, 600]);
    
    % 设置动画参数
    trail_length = 50; % 轨迹尾迹长度
    
    for i = 1:10:N
        clf;
        subplot(1, 2, 1);
        % XY平面轨迹
        start_idx = max(1, i-trail_length);
        plot(x_true(start_idx:i)*1000, y_true(start_idx:i)*1000, 'b-', 'LineWidth', 2); hold on;
        plot(x_est(start_idx:i)*1000, y_est(start_idx:i)*1000, 'r--', 'LineWidth', 2);
        plot(x_true(i)*1000, y_true(i)*1000, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
        plot(x_est(i)*1000, y_est(i)*1000, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
        xlabel('X位置 (mm)'); ylabel('Y位置 (mm)');
        title(sprintf('XY轨迹 - 时间: %.2f s', t(i)));
        legend('真实轨迹', '估计轨迹', 'Location', 'best');
        axis equal; grid on;
        
        % 设置坐标轴范围（自适应）
        all_x = [x_true, x_est];
        all_y = [y_true, y_est];
        x_range = [min(all_x)*1.1, max(all_x)*1.1] * 1000;
        y_range = [min(all_y)*1.1, max(all_y)*1.1] * 1000;
        xlim(x_range); ylim(y_range);
        
        subplot(1, 2, 2);
        % Rz轨迹
        plot(t(start_idx:i), rz_true(start_idx:i)*1000, 'b-', 'LineWidth', 2); hold on;
        plot(t(start_idx:i), rz_est(start_idx:i)*1000, 'r--', 'LineWidth', 2);
        plot(t(i), rz_true(i)*1000, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
        plot(t(i), rz_est(i)*1000, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
        xlabel('时间 (s)'); ylabel('Rz角度 (mrad)');
        title('Rz轨迹');
        legend('真实值', '估计值', 'Location', 'best');
        grid on;
        
        drawnow;
        pause(0.01);
    end
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