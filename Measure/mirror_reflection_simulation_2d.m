function h = cal_light_distance(para)
    h = zeros(3,1);
    h(1) = para.f(1) * para.g(1);
    h(2) = para.f(2) * para.g(2);
    h(3) = para.f(3) * para.g(3);
end

function main_func()
    max_errx = -1;
    max_erry = -1;
    max_errrz = -1;
    max_time = -1;
    
    tic;
    t_overhead = toc;

    t_slice = zeros(8,1);

    N = 10000;
    for iter = 1:N
        para.phi = [rand()*1e-6, rand()*1e-6, rand()*1e-6];
        para.psi = [rand()*1e-6, rand()*1e-6, rand()*1e-6];
        para.a = {[-0.2; 0.1], [-0.1; -0.2], [0.1; -0.2]};
        para.b = [0.1; 0.1; 0.1];
    
        x = 0 * 2 * (rand()-0.5); 
        y = 0 * rand(); 
        rz = 1e-4 * (rand()-0.5);
    
        para = update_para(x,y,rz,para);
        h = cal_light_distance(para);
        % [x_, y_, rz_] = cal_3dof(h, para, rz, 0);
        % 
        % x__ = h(1) - para.b(1);
        % y__ = mean([h(2)-para.b(2), h(3)-para.b(3)]);
        % rz__ = rz_;
        % 
        % % [xi, yi, rzi] = cal_3dof_iter(h, para, rz);
        tic;
        [xg, yg, rzg, t] = cal_3dof_grad(h, para, rz);
        elapsed_time_us = toc;
        t_slice = t_slice + t;

    
        % errx_ = x_ - x; erry_ = y_ - y; errrz_ = rz_ - rz;
        % errx__ = x__ - x;erry__ = y__ - y; errrz__ = rz__ - rz;
        % errxi = xi - x; erryi = yi - y; errrzi = rzi - rz;
        errxg = xg - x; erryg = yg - y; errrzg = rzg - rz;

        max_errx = max(abs(errxg), max_errx);
        max_erry = max(abs(erryg), max_erry);
        max_errrz = max(abs(errrzg), max_errrz);
        max_time = max(max_time, elapsed_time_us);
    end
    t_slice = (t_slice / N) * 1e6;
    t = (t) * 1e6;
    disp(t_slice);
    disp(t);
    disp(max_errx);
    disp(max_erry);
    disp(max_errrz);
    disp(max_time);
    return;
end

function [tx, ty, rz] = cal_3dof(h_measure, para, rz0, flag)
    Ly = para.a{2}(1) - para.a{1}(1);
    if flag == 1
        rz = rz0;
    else
        rz = atan((h_measure(3)-h_measure(2))/Ly);
        para = update_para(0, 0, rz, para);
    end

    U = zeros(3,2);
    v = zeros(3,1);
    for i = 1:3
        b = para.b(i);
        ax = para.a{i}(1);
        ay = para.a{i}(2);

        f = para.f(i);
        sin_zeta = para.sin_zeta(i);
        cos_zeta = para.cos_zeta(i);

        U(i,1) = f * sin_zeta; U(i,2) = -f * cos_zeta;
        v(i) = f * (b - ax * sin_zeta + ay * cos_zeta);
    end

    h_hat = h_measure - v;
    d = U \ h_hat;
    tx = d(1);
    ty = d(2);
end

function [tx, ty, rz, t] = cal_3dof_grad(h, para, rz0)
    t = zeros(8,1);
    tic;
    [x, y, rz0] = cal_3dof(h, para, rz0, 0);
    t(1) = toc;
    for i = 1:1
        tic;
        para = update_para(x, y, rz0, para);
        t(2) = toc;

        tic;
        grad = cal_gradient(para);
        t(3) = toc;

        tic;
        residual = cal_light_distance(para) - h;
        t(4) = toc;

        tic;
        dr = -50 * grad' * residual;
        t(5) = toc;

        tic;
        rz0 = rz0 + dr;
        t(6) = toc;

        tic;
        para = update_para(x, y, rz0, para);
        t(7) = toc;

        tic;
        [x, y, rz0] = cal_3dof(h, para, rz0, 1);
        t(8) = toc;
    end
    tx = x; ty = y; rz = rz0;
end

function grad = cal_gradient(para)
    grad = zeros(3,1);
    for i = 1:3
        beta = para.beta(i);
        f = para.f(i);

        sin_zeta = para.sin_zeta(i); 
        cos_zeta = para.cos_zeta(i);

        g = para.g(i);
        dx = para.dx(i);
        dy = para.dy(i);

        dfdr = -((2*beta^2+1)/(2*beta^2-1))^2 * sin(para.gamma(i));
        dgdr = dx * cos_zeta + dy * sin_zeta;

        grad(i) = dfdr * g + f * dgdr;
    end
end

function para = update_para(x, y, rz, para)
    para.gamma = zeros(3,1);
    para.zeta = zeros(3,1);
    para.beta = zeros(3,1);
    para.sin_zeta = zeros(3,1);
    para.cos_zeta = zeros(3,1);
    para.f = zeros(3,1);
    para.dx = zeros(3,1);
    para.dy = zeros(3,1);
    para.g = zeros(3,1);

    for i = 1:3
        para.gamma(i) = rz + para.psi(i) - para.phi(i);
        para.zeta(i) = rz + para.psi(i);
        if (i == 1)
            para.zeta(i) = para.zeta(i) - pi/2;
        end
        para.sin_zeta(i) = sin(para.zeta(i));
        para.cos_zeta(i) = cos(para.zeta(i));
        para.beta(i) = -cos(para.gamma(i));
        para.f(i) = para.beta(i) / (2*para.beta(i)^2 - 1);
        para.dx(i) = x - para.a{i}(1);
        para.dy(i) = y - para.a{i}(2);
        para.g(i) = para.b(i) + para.dx(i) * para.sin_zeta(i) - para.dy(i) * para.cos_zeta(i);
    end
    
end

main_func();