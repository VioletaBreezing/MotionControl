clear; clc; format long e

% --- 参数设置 ---
a = [0.0,   0.05, -0.05;
     0.05,  0.0,   0.0];
b = [0.1; 0.1; 0.1];
psi = [0; 1e-5; -1e-5];
phi = [0; 2e-5; -2e-5];

x_true = 0.012345678901234;
y_true = -0.008765432109876;
theta_true = 1.23456789e-3;

% 初值（带微小扰动）
x0 = x_true + 1e-10;
y0 = y_true - 2e-10;
theta0 = theta_true + 5e-7;

% --- 生成真实测量值 h ---
[h_true, ~] = sensor_model(x_true, y_true, theta_true, a, b, psi, phi);

% --- 单步 Newton-Raphson ---
x_est = x0;
y_est = y0;
theta_est = theta0;

% 计算当前估计下的 h 和 J
[h_est, J] = sensor_model(x_est, y_est, theta_est, a, b, psi, phi);

% 残差
r = h_true - h_est;

% 解 J * delta = r
delta = J \ r;

% 更新
x_new = x_est + delta(1);
y_new = y_est + delta(2);
theta_new = theta_est + delta(3);

% --- 误差计算 ---
err_x = abs(x_new - x_true);
err_y = abs(y_new - y_true);
err_theta = abs(theta_new - theta_true);

fprintf('Position error:\n');
fprintf('  |Δx| = %.3e m\n', err_x);
fprintf('  |Δy| = %.3e m\n', err_y);
fprintf('Angle error:\n');
fprintf('  |Δθ| = %.3e rad\n', err_theta);

if err_x < 1e-9 && err_y < 1e-9 && err_theta < 1e-9
    fprintf('\n✅ SUCCESS: All errors within 1 nm / 1 nrad!\n');
else
    fprintf('\n❌ FAILURE: Errors exceed tolerance.\n');
end