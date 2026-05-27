function h = cal_light_distance(x,y,rz,args)
    phi = args.phi;     % 干涉仪光线方向角
    psi = args.psi;     % 反射镜法线方向角
    a = args.a;         % 干涉仪中心位矢
    b = args.b;         % 反射镜到运动台坐标系原点距离

    gamma = zeros(3, 1);
    zeta = zeros(3, 1);
    dx = zeros(3, 1);
    dy = zeros(3, 1);

    for i = 1:3
        gamma(i) = rz + psi(i) - phi(i);
        zeta(i) = rz + psi(i);
        dx(i) = x - a{i}(1);
        dy(i) = y - a{i}(2);
    end

    zeta(1) = zeta(1) - pi/2;

    h = zeros(3, 1);
    for i = 1:3
        beta = -cos(gamma(i));
        f = beta / (2*beta^2 - 1);
        g = b(i) + dx(i) * sin(zeta(i)) - dy(i) * cos(zeta(i));
        h(i) = f * g;
    end
end

function main_func()
    para.phi = [1e-6, 1e-6, -1e-6];
    para.psi = [1e-6, -1e-6, 1e-6];
    para.a = {[-0.2; 0.1], [-0.1; -0.2], [0.1; -0.2]};
    para.b = [0.1; 0.1; 0.1];

    x = 0; y = 0; rz = 0;
    h = cal_light_distance(x,y,rz,para);

    h
end

main_func();