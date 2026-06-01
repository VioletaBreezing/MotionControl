clear; clc; close all;

% --- Nominal (calibrated) model parameters used in ESTIMATION ---
a_nom = [0.0,   0.05, -0.05;
         0.05,  0.0,   0.0];
b_nom = [0.1; 0.1; 0.1];
psi_nom = [0; 1e-5; -1e-5];      % What estimator believes
phi_nom = [0; 2e-5; -2e-5];

% --- Actual (true) physical parameters (10% mismatch) ---
psi_actual = 1.1 * psi_nom;      % +10% error in psi
phi_actual = 0.9 * phi_nom;      % -10% error in phi
% Note: a and b are geometry, assumed perfect (or calibrated out)

% --- Trajectory & noise ---
fs = 1 / 2e-4;          % 5000 Hz
dt = 1/fs;
t = 0:dt:10;
N = length(t);

x_true = 2e-3 * sin(2*pi*t);
y_true = 1e-1 * sin(2*pi*t);      % ±100 mm!
theta_true = 2e-3 * sin(2*pi*t);

noise_rms = 0.31e-9;    % 0.31 nm RMS
rng(42); % reproducibility

% Preallocate
x_est = zeros(size(t));
y_est = zeros(size(t));
theta_est = zeros(size(t));
err_x = zeros(size(t));
err_y = zeros(size(t));
err_theta = zeros(size(t));

% Initial guess
x_est(1) = 0.0;
y_est(1) = 0.0;
theta_est(1) = 0.0;

fprintf('Running estimation with:\n');
fprintf('  - %.2f nm RMS sensor noise\n', noise_rms*1e9);
fprintf('  - 10%% mismatch in psi/phi (estimator uses nominal model)\n');

for k = 1:N
    % --- TRUE WORLD: generate measurement with ACTUAL parameters ---
    [h_true, ~] = sensor_model(x_true(k), y_true(k), theta_true(k), ...
                               a_nom, b_nom, psi_actual, phi_actual);
    h_meas = h_true + noise_rms * randn(3, 1); % add noise
    
    % --- ESTIMATOR: uses NOMINAL model (mismatched!) ---
    if k == 1
        x0 = x_est(1); y0 = y_est(1); theta0 = theta_est(1);
    else
        x0 = x_est(k-1); y0 = y_est(k-1); theta0 = theta_est(k-1);
    end
    
    x_cur = x0; y_cur = y0; theta_cur = theta0;
    for iter = 1:2
        % ⚠️ Here we use NOMINAL model to compute h_est and J!
        [h_est, J] = sensor_model(x_cur, y_cur, theta_cur, ...
                                  a_nom, b_nom, psi_nom, phi_nom);
        r = h_meas - h_est;
        delta = J \ r;
        x_cur = x_cur + delta(1);
        y_cur = y_cur + delta(2);
        theta_cur = theta_cur + delta(3);
    end
    
    x_est(k) = x_cur;
    y_est(k) = y_cur;
    theta_est(k) = theta_cur;
    
    err_x(k) = x_cur - x_true(k);
    err_y(k) = y_cur - y_true(k);
    err_theta(k) = theta_cur - theta_true(k);
end

% --- Analysis: separate BIAS and NOISE ---
bias_x = mean(err_x);
bias_y = mean(err_y);
bias_theta = mean(err_theta);

rms_noise_x = sqrt(mean((err_x - bias_x).^2));
rms_noise_y = sqrt(mean((err_y - bias_y).^2));
rms_noise_theta = sqrt(mean((err_theta - bias_theta).^2));

fprintf('\n✅ FINAL RESULT (10%% model mismatch + 0.31 nm noise):\n');
fprintf('BIAS   : Δx = %+.3f nm, Δy = %+.3f nm, Δθ = %+.3f nrad\n', ...
        bias_x*1e9, bias_y*1e9, bias_theta*1e9);
fprintf('NOISE  : σ_x = %.3f nm, σ_y = %.3f nm, σ_θ = %.3f nrad\n', ...
        rms_noise_x*1e9, rms_noise_y*1e9, rms_noise_theta*1e9);

% --- Plot ---
figure;
subplot(3,1,1);
plot(t, err_x*1e9, 'r'); ylabel('\Delta x (nm)'); grid on;
title('Estimation Error with Model Mismatch + Noise');
hold on; plot([t(1), t(end)], [bias_x*1e9, bias_x*1e9], 'k--');
subplot(3,1,2);
plot(t, err_y*1e9, 'g'); ylabel('\Delta y (nm)'); grid on;
hold on; plot([t(1), t(end)], [bias_y*1e9, bias_y*1e9], 'k--');
subplot(3,1,3);
plot(t, err_theta*1e9, 'b'); ylabel('\Delta \theta (nrad)'); xlabel('Time (s)'); grid on;
hold on; plot([t(1), t(end)], [bias_theta*1e9, bias_theta*1e9], 'k--');
legend('Error','Bias');

% --- Assessment ---
fprintf('\n📊 Interpretation:\n');
if abs(bias_x)*1e9 < 0.5 && abs(bias_y)*1e9 < 0.5 && abs(bias_theta)*1e9 < 0.5
    fprintf('🎉 BIAS is negligible (< 0.5 nm/nrad) — system robust to 10%% mismatch!\n');
else
    fprintf('⚠️  Significant BIAS detected — consider online recalibration.\n');
end