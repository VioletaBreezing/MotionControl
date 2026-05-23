function phi = SO3_to_so3_my(R)
    % 将 SO(3) 中的旋转矩阵 R 转换为 so(3) 中的旋转向量 phi
    % 输入: R (3x3 旋转矩阵)
    % 输出: phi (3x1 向量，单位：弧度)

    % 确保 R 是有效的旋转矩阵（可选）
    assert(all(size(R) == [3,3]), 'R must be 3x3');
    
    % 计算旋转角 theta
    cos_theta = (trace(R) - 1) / 2;
    cos_theta = max(-1.0, min(1.0, cos_theta)); % 防止数值误差导致 acos 越界
    theta = acos(cos_theta);

    if abs(theta) < 1e-10
        % 情况1: 无旋转 → phi = [0;0;0]
        phi = zeros(3,1);
    elseif theta > pi - 1e-6
        % 情况2: 接近180度旋转（数值不稳定，需特殊处理）
        % 利用 R + I 的最大对角元确定主轴
        A = (R + eye(3)) / 2;
        [~, idx] = max(diag(A));
        n = A(:, idx) / norm(A(:, idx));
        phi = pi * n;
    else
        % 情况3: 一般情况
        sin_theta = sin(theta);
        nx = (R(3,2) - R(2,3)) / (2*sin_theta);
        ny = (R(1,3) - R(3,1)) / (2*sin_theta);
        nz = (R(2,1) - R(1,2)) / (2*sin_theta);
        n = [nx; ny; nz];
        phi = theta * n;
    end
end