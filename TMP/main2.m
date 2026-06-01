clear; clc; format long e

% --- 固定参数 ---
a = [0.0,   0.05, -0.05;
     0.05,  0.0,   0.0];
b = [0.1; 0.1; 0.1];
psi = [0; 1e-5; -1e-5];
phi = [0; 2e-5; -2e-5];

x_true = 0.012345678901234;
y_true = -0.008765432109876;
theta_true = 1.23456789e-3;

% 生成真实测量值
[h_true, ~] = sensor_model(x_true, y_true, theta_true, a, b, psi, phi);

% --- 扰动设置 ---
dx_list = [1e-6, 1e-5, 1e-4];      % 1μm, 10μm, 100μm
dtheta_list = [1e-3, 5e-3, 1e-2];  % 1mrad, 5mrad, 10mrad

fprintf('Testing robustness to initial perturbations...\n');
fprintf('---------------------------------------------\n');

for dx = dx_list
    for dtheta = dtheta_list
        % 初值（最坏情况：x,y 同向最大扰动，theta 也最大）
        x0 = x_true + dx;
        y0 = y_true + dx;  % same magnitude
        theta0 = theta_true + dtheta;

        % --- Newton-Raphson (最多2步) ---
        x_est = x0; y_est = y0; theta_est = theta0;
        converged = false;
        
        for iter = 1:2
            [h_est, J] = sensor_model(x_est, y_est, theta_est, a, b, psi, phi);
            r = h_true - h_est;
            
            % Check residual norm
            res_norm = norm(r);
            if res_norm < 1e-15
                converged = true;
                break;
            end
            
            delta = J \ r;
            x_est = x_est + delta(1);
            y_est = y_est + delta(2);
            theta_est = theta_est + delta(3);
        end

        err_x = abs(x_est - x_true);
        err_y = abs(y_est - y_true);
        err_theta = abs(theta_est - theta_true);
        total_err = max([err_x, err_y, err_theta]);

        status = '✅';
        if total_err > 1e-9
            status = '⚠️ ';
        end
        if any(isnan([x_est, y_est, theta_est]))
            status = '❌';
        end

        fprintf('%s Δx₀=%.0fμm, Δθ₀=%.1fmrad → err=%.1e m/rad (iters=%d)\n', ...
                status, dx*1e6, dtheta*1e3, total_err, iter);
    end
end

fprintf('\n✅ = within 1nm/1nrad, ⚠️ = exceeds but finite, ❌ = NaN/diverge\n');