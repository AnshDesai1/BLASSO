% Sliding Frank-Wolfe Algorithm (Algorithm 2) in MATLAB
% Solves: min_m 1/2 * ||Phi m - y||^2 + lambda*||m||_{TV} for spike deconvolution
% https://arxiv.org/pdf/1811.06416
function [theta_est, a_est] = sliding_frank_wolfe(y_obs, y_grid, d, lambda, max_iter, tol)
    if nargin < 6, tol = 1e-6; end
    if nargin < 5, max_iter = 50; end

    theta_est = [];
    a_est = [];
    for k = 1:max_iter
        if isempty(theta_est)
            r = y_obs;  % m^(0) = 0
        else
            A_cur = exp(2i * pi * d * (y_grid * sin(theta_est')));
            r = y_obs - A_cur * a_est;
        end
        % Compute dual certificate eta = 1/lambda * Phi^*(y-Phi*m)
        eta_fn = @(theta) (1/lambda) * exp(-2i * pi * d * sin(theta) * y_grid') * r;
        obj_fn = @(theta) -abs(eta_fn(theta));
        theta_star = fminbnd(obj_fn, -pi, pi);
        wrap_diff = abs(mod(theta_est - theta_star + pi, 2*pi) - pi);
        if any(wrap_diff < 1e-2)
            continue;
        end

        if abs(eta_fn(theta_star)) <= 1- 1e-2
            break;
        end

        % Augment support and solve LASSO
        theta_est_half = [theta_est; theta_star];
        A_half = exp(2i * pi * d * (y_grid * sin(theta_est_half')));
        a_est_half = solve_lasso_complex(A_half, y_obs, lambda);

        N = length(theta_est_half);
        z0 = [real(a_est_half); imag(a_est_half); theta_est_half];
        noise_level = 0.01 / k;
        z0 = z0 + noise_level * randn(size(z0));

        % Lower and upper bounds (only constrain theta, not amplitudes)
        lb = [-Inf(2*N,1); -pi/2 * ones(N,1)];
        ub = [ Inf(2*N,1);  pi/2 * ones(N,1)];
        
        opts = optimoptions('fmincon', ...
            'Algorithm', 'interior-point', ...
            'Display', 'none', ...
            'OptimalityTolerance', 1e-10, ...
            'StepTolerance', 1e-10, ...
            'FunctionTolerance', 1e-10, ...
            'MaxIterations', 1000, ...
            'MaxFunctionEvaluations', 5000);
        
        cost_fn = @(z) sliding_objective(z, y_obs, y_grid, d, lambda);
        [z_opt, fval, exitflag] = fmincon(cost_fn, z0, [], [], [], [], lb, ub, [], opts);

        a_est = z_opt(1:N) + 1i * z_opt(N+1:2*N);
        theta_est = z_opt(2*N+1:end);

        mask = abs(a_est) > 1e-3;
        a_est = a_est(mask);
        theta_est = theta_est(mask);

        [theta_est, ia] = uniquetol(theta_est, 1e-3);
        a_est = a_est(ia);

        a_opt = a_est;
        A_opt = exp(2i * pi * d * (y_grid * sin(theta_est')));
        r = y_obs - A_opt * a_opt;
        
        loss = 0.5 * norm(r)^2 + lambda * sum(abs(a_opt));
        fprintf('Iter %2d: Loss = %.6f\n', k, loss);
    end
    % Pruning
    if exist('a_opt', 'var') && ~isempty(a_opt) % Checking if loop did not quit instantly
        theta_range = [-pi/2, pi/2]; % Restrict theta range
        mask = abs(a_opt) > 0.001 & theta_est >= theta_range(1) & theta_est <= theta_range(2); % Pruning bad guesses
        theta_est = theta_est(mask);
        a_est = a_opt(mask);
    else
        theta_est = [];
        a_est = [];
        warning('No spikes selected. Empty estimate returned.');
    end
end

function val = sliding_objective(z, y_obs, y_grid, d, lambda)
    N = length(z) / 3;
    a_real = z(1:N);
    a_imag = z(N+1:2*N);
    theta = z(2*N+1:end);
    a = a_real + 1i * a_imag;
    A = exp(2i * pi * d * (y_grid * sin(theta')));
    res = A * a - y_obs;
    val = 0.5 * norm(res)^2 + lambda * sum(abs(a));
end

function a = solve_lasso_complex(A, X, lambda)
    % Realify complex LASSO: A a = X becomes [Re A, -Im A; Im A, Re A][Re a; Im a] = [Re X; Im X]
    A_real = [real(A), -imag(A); imag(A), real(A)];
    X_real = [real(X); imag(X)];
    N = size(A, 2);
    cvx_begin quiet
        variable t_real(2*N)
        minimize(0.5 * sum_square(A_real * t_real - X_real) + lambda * norm(t_real, 1))
    cvx_end
    a = t_real(1:N) + 1i * t_real(N+1:end);
end