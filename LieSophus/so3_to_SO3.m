function R = so3_to_SO3(phi)
    theta = norm(phi);
    if theta < 1e-8
        R = eye(3);
    else
        k = phi / theta;
        K = skew(k);
        % R = eye(3) + sin(theta)*K + (1-cos(theta))*K*K;
        R = cos(theta)*eye(3) + (1-cos(theta))*(k*k') + sin(theta)*K;
    end
end
