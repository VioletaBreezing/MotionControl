function J = J_left(phi)
    theta = norm(phi);
    if theta < 1e-10
        J = eye(3);
    else
        a = phi / theta;
        J = sin(theta)/theta * eye(3) + (1 - sin(theta)/theta) * (a * a') + ...
            (1 - cos(theta))/theta * skew(a);
    end
end
