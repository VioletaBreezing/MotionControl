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

% --- Sensor noise (typical high-end interferometer) ---
noise_rms = 0.31e-9;    % 0.31 nanometer RMS
rng(42); % for reproducibility

% --- Initial guess ---
x_est(1) = 0.0;
y_est(1) = 0.0;
theta_est(1) = 0.0;

fprintf('Running closed-loop estimation with %.2f nm RMS noise...\n', noise_rms*1e9);

for k = 1:N
    % Generate true measurement
    [h_true, ~] = sensor_model(x_true(k), y_true(k), theta_true(k), a, b, psi, phi);
    
    % Add realistic sensor noise
    h_meas = h_true + noise_rms * randn(3, 1);
    
    % Initial guess: use previous estimate
    if k == 1
        x0 = x_est(1);
        y0 = y_est(1);
        theta0 = theta_est(1);
    else
        x0 = x_est(k-1);
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
    
    % Store estimate
    x_est(k) = x_cur;
    y_est(k) = y_cur;
    theta_est(k) = theta_cur;
    
    % Compute error
    err_x(k) = x_cur - x_true(k);  % keep signed for RMS
    err_y(k) = y_cur - y_true(k);
    err_theta(k) = theta_cur - theta_true(k);
end

% --- Final analysis: use RMS error ---
rms_err_x = sqrt(mean(err_x.^2));
rms_err_y = sqrt(mean(err_y.^2));
rms_err_theta = sqrt(mean(err_theta.^2));

baseline = 0.05; % meters (from your 'a' matrix: max baseline ~5 cm)
equiv_noise_theta = noise_rms / baseline; % rad

fprintf('\n✅ FINAL RESULT (with %.2f nm RMS sensor noise):\n', noise_rms*1e9);
fprintf('RMS |Δx| = %.3f nm\n', rms_err_x*1e9);
fprintf('RMS |Δy| = %.3f nm\n', rms_err_y*1e9);
fprintf('Equivalent angle noise floor: %.2f nrad\n', equiv_noise_theta*1e9);
fprintf('Achieved RMS Δθ             : %.2f nrad\n', rms_err_theta*1e9);

if rms_err_x < 2*noise_rms && ...
   rms_err_y < 2*noise_rms && ...
   rms_err_theta < 2*equiv_noise_theta
    fprintf('\n🎉 SUCCESS: All errors near theoretical limit!\n');
else
    fprintf('\n⚠️  Check model or noise assumptions.\n');
end
% --- Plot results ---
figure;
subplot(3,1,1);
plot(t, x_true*1e9, 'k', t, x_est*1e9, '--r'); ylabel('x (nm)'); legend('True','Est'); grid on;
title(sprintf('Pose Estimation with %.2f nm RMS Noise', noise_rms*1e9));
subplot(3,1,2);
plot(t, y_true*1e6, 'k', t, y_est*1e6, '--r'); ylabel('y (μm)'); legend('True','Est'); grid on;
subplot(3,1,3);
plot(t, theta_true*1e6, 'k', t, theta_est*1e6, '--r'); ylabel('\theta (μrad)'); xlabel('Time (s)'); legend('True','Est'); grid on;

figure;
plot(t, err_x*1e9, 'r', t, err_y*1e9, 'g', t, err_theta*1e6, 'b');
xlabel('Time (s)'); ylabel('Error'); 
legend('\Delta x (nm)','\Delta y (nm)','\Delta \theta (\mu rad)'); grid on;
title('Estimation Error Over Time (with Noise)');