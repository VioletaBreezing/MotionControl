function [x, y, rz] = cal_eddy_3dof(hx, hy1, hy2, param)

    ax = param.ax; ay1 = param.ay1; ay2 = param.ay2;
    bx = param.bx; by = param.by;
    ay = ay1 + ay2;

    x = hx;
    y = ay2/ay * hy1 + ay1/ay * hy2;
    rz = atan((hy2 - hy1) / ay);

    for iter = 1:20
        armx  = -y + ax  - bx * tan(rz/2);
        army1 =  x + ay1 + by * tan(rz/2);
        army2 = -x + ay2 - by * tan(rz/2);
    
        abbey_errx  = armx  * tan(rz);
        abbey_erry1 = army1 * tan(rz);
        abbey_erry2 = army2 * tan(rz);
    
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

% param.ax = 0.09;
% param.ay1 = 0.107;
% param.ay2 = 0.107;
% param.bx = 0.1975;
% param.by = 0.116;
% 
% N = 10000;
% hx_v  = 1e-3 * randn(N, 1);
% hy1_v = 1e-3 * randn(N, 1);
% hy2_v = 1e-3 * randn(N, 1);
% pos = zeros(3, N);
% 
% for i = 1:N
%     [x,y,rz] = cal_eddy_3dof(hx_v(i), hy1_v(i), hy2_v(i), param);
%     pos(1, i) = x;
%     pos(2, i) = y;
%     pos(3, i) = rz;
% end
% 
% figure;
% for i = 1:3
%     subplot(3,1,i);
%     plot(pos(i, :));
%     grid on;
% end

function y = my_tan(x)
    y = x + x^3/3;
end

function y = my_atan(x)
    y = x - x^3/3;
end