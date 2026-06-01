%% 传感器测量模型参数标定仿真（Levenberg-Marquardt）
clear; close all; clc;

% -----------------------------
% 1. 真实参数（用于生成仿真数据）
% -----------------------------
phi_true = [1.2e-3, -0.8e-3, 0.5e-3];     % rad
psi_true = [-0.9e-3, 1.1e-3, -0.6e-3];    % rad

a_true = [ -0.2,  0.1;   % a0x, a0y
           -0.1, -0.2;   % a1x, a1y
            0.1, -0.2 ]; % a2x, a2y   (单位: m)

b_true = [0.1; 0.1; 0.1];         % m

% 打包真实参数向量 p = [phi0,phi1,phi2, psi0,psi1,psi2, a0x,a0y,a1x,a1y,a2x,a2y, b0,b1,b2]
p_true = [phi_true, psi_true, a_true(:)', b_true'];

% -----------------------------
% 2. 生成仿真输入数据 (x, y, theta)
% -----------------------------
N = 50/2e-2;  % 数据点数

% x: -1mm ~ +1mm
x_data = (-1e-3 + 2e-3 * rand(N,1));  
% y: -0.5cm ~ 23cm → -0.005m ~ 0.23m
y_data = (-0.005 + (0.23 + 0.005) * rand(N,1));
% theta: -2mrad ~ +2mrad
theta_data = (-2e-3 + 4e-3 * rand(N,1));

% 引入少量刻意变化以激发角度参数（可选增强）
theta_data(1:10) = linspace(-2e-3, 2e-3, 10)';

x_data = (-0.25 + 0.5 * rand(N,1));      % -250mm ~ +250mm （覆盖 a_ix 范围）
y_data = (-0.1 + 0.4 * rand(N,1));       % -100mm ~ +300mm
theta_data = linspace(-5e-3, 5e-3, N)';  % ±5 mrad

% x: -1mm ~ +1mm
t = (0:N-1) * 2e-2;
x_data = 2e-3 * sin(2*pi*2*t);  
% y: -0.5cm ~ 23cm → -0.005m ~ 0.23m
y_data = 1e-1 * sin(2*pi*t);
% theta: -2mrad ~ +2mrad
theta_data = 2e-3 * sin(2*pi*4*t);

% -----------------------------
% 3. 生成无噪声和带噪声的测量值 h
% -----------------------------
h_clean = zeros(N,3);
h_meas  = zeros(N,3);
noise_rms = 0.31e-9/6;

for k = 1:N
    x = x_data(k); y = y_data(k); th = theta_data(k);
    h_clean(k,:) = sensor_model(x, y, th, p_true)';
end

% 添加高斯噪声
h_meas = h_clean + noise_rms * randn(N,3);

% -----------------------------
% 4. 初始猜测（基于先验知识）
% -----------------------------
phi0_guess = [0, 0, 0];          % 1x3
psi0_guess = [0, 0, 0];          % 1x3
a0_guess = [ -0.22,  0.18;
             -0.11, -0.26;
              0.15, -0.24 ];       % 3x2
% b0_guess = mean(h_meas, 1);      % 1x3 （关键：不要转置！）
% 初始猜测改进
b0_guess = -mean(h_meas, 1);  % 注意负号！

p0 = [phi0_guess, psi0_guess, a0_guess(:)', b0_guess];  % 全是 1x?，可拼接

% -----------------------------
% 5. 设置优化选项（Levenberg-Marquardt）
% -----------------------------
opts = optimoptions('lsqnonlin', ...
    'Algorithm', 'levenberg-marquardt', ...
    'Display', 'iter', ...
    'MaxIterations', 500, ...
    'FunctionTolerance', 1e-12, ...
    'StepTolerance', 1e-12);

% 将数据打包传入残差函数
data = struct('x', x_data, 'y', y_data, 'theta', theta_data, 'h_meas', h_meas);

% -----------------------------
% 6. 执行优化
% -----------------------------
[p_est, resnorm, residual, exitflag, output] = lsqnonlin(@(p) residual_func(p, data), p0, [], [], opts);

% -----------------------------
% 7. 结果分析
% -----------------------------
fprintf('\n=== 标定结果 ===\n');
param_names = {'phi0','phi1','phi2','psi0','psi1','psi2',...
               'a0x','a0y','a1x','a1y','a2x','a2y','b0','b1','b2'};
for i = 1:length(p_true)
    fprintf('%5s: 真实值 = %10.4e, 估计值 = %10.4e, 误差 = %10.4e\n', ...
        param_names{i}, p_true(i), p_est(i), p_est(i)-p_true(i));
end

fprintf('\n残差平方和: %.3e\n', resnorm);
fprintf('退出标志: %d (%s)\n', exitflag, output.message);

% 可视化：预测 vs 测量
h_pred = zeros(N,3);
for k = 1:N
    h_pred(k,:) = sensor_model(x_data(k), y_data(k), theta_data(k), p_est)';
end

figure;
for i = 1:3
    subplot(3,1,i);
    plot(1:N, h_meas(:,i), 'r.', 'MarkerSize', 8); hold on;
    plot(1:N, h_pred(:,i), 'b-', 'LineWidth', 1.2);
    plot(1:N, h_clean(:,i), 'k--', 'LineWidth', 1);
    legend('测量值','模型预测','真值');
    ylabel(['h_' num2str(i-1)]);
    grid on;
end
xlabel('样本索引');
sgtitle('传感器读数：真值、测量值与模型预测对比');

% -----------------------------
% 8. 函数定义
% -----------------------------

function h = sensor_model(x, y, theta, p)
    % 解包参数
    phi = p(1:3);
    psi = p(4:6);
    a = reshape(p(7:12), [3,2]);  % [a0x,a0y; a1x,a1y; a2x,a2y]
    b = p(13:15);

    h = zeros(3,1);
    for i = 1:3
        gamma = theta + psi(i) - phi(i);
        f = -cos(gamma) / cos(2*gamma);
        
        if i == 1  % i=0 in math, but MATLAB index starts at 1
            zeta = theta + psi(1) - pi/2;
        else
            zeta = theta + psi(i);
        end
        
        g = b(i) + (x - a(i,1))*sin(zeta) - (y - a(i,2))*cos(zeta);
        h(i) = f * g;
    end
end

function res = residual_func(p, data)
    % 残差函数：返回 (N*3)x1 向量
    x = data.x; y = data.y; theta = data.theta; h_meas = data.h_meas;
    N = length(x);
    res_vec = zeros(N*3, 1);
    
    for k = 1:N
        h_pred = sensor_model(x(k), y(k), theta(k), p);
        res_vec((k-1)*3 + (1:3)) = h_pred - h_meas(k,:)';
    end
    res = res_vec;
end