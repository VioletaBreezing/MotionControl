function T = se3_to_SE3(xi)
    if nargin ~= 1
        error('需要一个输入参数');
    end
    
    % 情况1：输入是 6x1 向量
    if isvector(xi) && length(xi) == 6
        rho = xi(1:3);
        phi = xi(4:6);
        xi_hat = [skew(phi), rho; 0 0 0 0];  % 构造 4x4 se(3) 矩阵
    elseif size(xi,1) == 4 && size(xi,2) == 4
        xi_hat = xi;
        % 从 xi_hat 提取 phi 和 rho
        phi = [xi(3,2); xi(1,3); xi(2,1)];  % 因为 skew(phi) = [0 -w3 w2; w3 0 -w1; -w2 w1 0]
        rho = xi(1:3,4);
    else
        error('输入必须是6x1向量或4x4矩阵');
    end
    
    % 计算旋转角度
    theta = norm(phi);
    
    % 初始化 R 和 rho
    if theta < 1e-8
        % 小角度近似：R ≈ I, t ≈ [1;1;1]
        R = eye(3);
        t = rho;
    else
        R = so3_to_SO3(phi);
        
        % 计算 t（用于平移部分）
        J_ = J_left(phi);
        t = J_ * rho;
    end
    
    % 构造 SE(3) 矩阵
    T = zeros(4);
    T(1:3,1:3) = R;
    T(1:3,4)   = t;
    T(4,4) = 1;
end