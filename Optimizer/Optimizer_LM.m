function [x_hat, resnorm, converged, iter] = Optimizer_LM(...
        ResFcn, Rargs, JacobiFcn, Jargs, ...
        x0, ...
        TolX, TolR, TolLamda, MaxIter, ...
        debug)

    lambda = 1e-3;
    nu = 2;
    x = x0;
    converged = false;

    for iter = 1:MaxIter
        r = ResFcn(x, Rargs);
        J = JacobiFcn(x, Jargs);
        cost = 0.5 * sum(r.^2);

        A = J' * J + lambda * eye(6);
        b = -J' * r;
        dx = A \ b;
        
        if debug
            fprintf('[%s] Iter %2d: cost = %.3e\n', mfilename, iter, cost);
        end

        if norm(r) < TolR
            if debug
                fprintf('[%s]: Converged by residual norm.\n', mfilename);
            end
            converged = true;
            break;
        end

        if norm(dx) < TolX
            if norm(r) > TolR
                fprintf('[%s] LM failed: Residual(%.3e) norm cannot converge.', mfilename, norm(r));
            end
            if debug
                fprintf('[%s]: Ended iteration by X change.\n', mfilename);
            end
            break;
        end

        x_new = x + dx;
        r_new = ResFcn(x_new, Rargs);
        cost_new = 0.5 * sum(r_new.^2);

        rho = (cost - cost_new) / (0.5 * dx' * (lambda * dx - J' * r));

        % 接受或拒绝更新
        if rho > 0
            x = x_new;
            lambda = lambda * max(1/3, 1 - (2*rho - 1)^3);
            nu = 2;
        else
            lambda = lambda * nu;
            nu = 2 * nu;
        end
        
        if lambda > TolLamda
            error('[%s] failed: lambda = %.4e too large.', mfilename, lambda);
        end
    end
    x_hat = x;
    resnorm = r;
end
