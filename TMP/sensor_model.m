function [h, J] = sensor_model(x, y, theta, a, b, psi, phi)
    % 输入：标量 x, y, theta；a(2,3), b(3,1), psi(3,1), phi(3,1)
    % 输出：h(3,1), J(3,3) = [dh/dx, dh/dy, dh/dtheta]

    h = zeros(3,1);
    J = zeros(3,3);

    for i = 1:3
        gamma_i = theta + psi(i) + phi(i);
        zeta_i = theta + psi(i);
        if i == 1
            zeta_i = zeta_i - pi/2;  % i=0 in your notation (MATLAB index 1)
        end

        cos_gamma = cos(gamma_i);
        cos_2gamma = cos(2*gamma_i);
        f_i = -cos_gamma / cos_2gamma;

        sin_zeta = sin(zeta_i);
        cos_zeta = cos(zeta_i);

        dx = x - a(1,i);
        dy = y - a(2,i);

        g_i = b(i) + dx * sin_zeta - dy * cos_zeta;
        h(i) = f_i * g_i;

        % --- Jacobian ---
        % df_i / dtheta
        df_dtheta = (sin(gamma_i)*cos(2*gamma_i) - 2*cos(gamma_i)*sin(2*gamma_i)) / (cos_2gamma^2);
        
        % d(sin_zeta)/dtheta = cos(zeta_i)
        % d(cos_zeta)/dtheta = -sin(zeta_i)
        % (无论是否减 pi/2，导数形式不变)
        dsin_dtheta = cos_zeta;
        dcos_dtheta = -sin_zeta;

        dg_dx = sin_zeta;
        dg_dy = -cos_zeta;
        dg_dtheta = dx * dsin_dtheta - dy * dcos_dtheta;

        dh_dx = f_i * dg_dx;
        dh_dy = f_i * dg_dy;
        dh_dtheta = df_dtheta * g_i + f_i * dg_dtheta;

        J(i, :) = [dh_dx, dh_dy, dh_dtheta];
    end
end