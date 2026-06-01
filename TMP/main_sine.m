clear; clc; close all;

% --- 系统参数 ---
a = [0.0,   0.05, -0.05;
     0.05,  0.0,   0.0];
b = [0.1; 0.1; 0.1];
psi = [0; 1e-5; -1e-5];
phi = [0; 2e-5; -2e-5];

% --- 轨迹设置 ---
fs = 1 / 2e-4;          % 5000 Hz
dt = 1/fs;
t = 0:dt:10;            % 10 seconds
N = length(t);

% True trajectory (1 Hz sine)
x_true = 2e-3 * sin(2*pi*t);
y_true = 1e-1 * sin(2*pi*t);      % ±100 mm!
theta_true = 2e-3 * sin(2*pi*t);

% Preallocate
x_est = zeros(size(t));
y_est = zeros(size(t));
theta_est = zeros(size(t));
err_x = zeros(size(t));
err_y = zeros(size(t));
err_theta = zeros(size(t));

% --- Initial guess at t=0 (realistic: e.g., from calibration or power-on state) ---
x_est(1) = 0.0;           % Could also be small offset, e.g., 1e-3
y_est(1) = 0.0;
theta_est(1) = 0.0;

fprintf('Running closed-loop estimation (corrected)...\n');

for k = 1:N
    % Get measurement from TRUE state (sensor sees reality)
    [h_meas, ~] = sensor_model(x_true(k), y_true(k), theta_true(k), a, b, psi, phi);

    % Optional: add realistic sensor noise (e.g., 0.31 nm RMS)
    h_meas = h_meas + 0.31e-9 * randn(3,1);
    
    % --- Set initial guess for Newton-Raphson ---
    if k == 1
        x0 = x_est(1);      % initial guess
        y0 = y_est(1);
        theta0 = theta_est(1);
    else
        x0 = x_est(k-1);    % use PREVIOUS ESTIMATE as initial guess
        y0 = y_est(k-1);
        theta0 = theta_est(k-1);
    end
    
    % --- 2-step Newton-Raphson ---
    x_cur = x0; y_cur = y0; theta_cur = theta0;
    for iter = 1:2
        [h_est, J] = sensor_model(x_cur, y_cur, theta_cur, a, b, psi, phi);
        r = h_meas - h_est;
        delta = J \ r;
        x_cur = x_cur + delta(1);
        y_cur = y_cur + delta(2);
        theta_cur = theta_cur + delta(3);
    end
    
    % Store current estimate
    x_est(k) = x_cur;
    y_est(k) = y_cur;
    theta_est(k) = theta_cur;
    
    % Compute error
    err_x(k) = abs(x_cur - x_true(k));
    err_y(k) = abs(y_cur - y_true(k));
    err_theta(k) = abs(theta_cur - theta_true(k));
    
    if mod(k, 5000) == 0 || k == 1
        fprintf('  Step %d/%d: |Δx|=%.2e m, |Δy|=%.2e m, |Δθ|=%.2e rad\n', ...
                k, N, err_x(k), err_y(k), err_theta(k));
    end
end

% --- Final analysis ---
max_err_x = max(err_x);
max_err_y = max(err_y);
max_err_theta = max(err_theta);

fprintf('\n✅ FINAL RESULT (Closed-loop, correct initial guess propagation):\n');
fprintf('Max |Δx| = %.3e m  (%.3f nm)\n', max_err_x, max_err_x*1e9);
fprintf('Max |Δy| = %.3e m  (%.3f nm)\n', max_err_y, max_err_y*1e9);
fprintf('Max |Δθ| = %.3e rad (%.3f nrad)\n', max_err_theta, max_err_theta*1e9);

if max_err_x < 1e-9 && max_err_y < 1e-9 && max_err_theta < 1e-9
    fprintf('\n🎉 SUCCESS: All errors within 1 nm / 1 nrad in realistic closed-loop!\n');
else
    fprintf('\n⚠️  WARNING: Errors exceed tolerance.\n');
end

% --- Plot ---
figure;
subplot(3,1,1);
plot(t, x_true*1e3, 'k', t, x_est*1e3, '--r'); ylabel('x (mm)'); legend('True','Est'); grid on;
title('Closed-loop Pose Estimation (Corrected)');
subplot(3,1,2);
plot(t, y_true*1e3, 'k', t, y_est*1e3, '--r'); ylabel('y (mm)'); legend('True','Est'); grid on;
subplot(3,1,3);
plot(t, theta_true*1e3, 'k', t, theta_est*1e3, '--r'); ylabel('\theta (mrad)'); xlabel('Time (s)'); legend('True','Est'); grid on;

figure;
semilogy(t, err_x, 'r', t, err_y, 'g', t, err_theta, 'b');
xlabel('Time (s)'); ylabel('Absolute Error (log scale)'); 
legend('|Δx|','|Δy|','|Δθ|'); grid on;
title('Estimation Error Over Time');