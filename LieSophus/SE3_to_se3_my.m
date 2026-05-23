function xi = SE3_to_se3_my(T)
    assert(all(size(T) == [4,4]), 'T must be 4x4');

    R = T(1:3,1:3);
    t = T(1:3, 4);

    phi = SO3_to_so3_my(R);
    J = J_left(phi);
    rho = J\t;
    xi = [rho; phi];
end
