function [phi_new, residuals, num_iters] = column_sweep(phi_old, source, deltax_sq, deltay_sq, tol, max_iter)
    if nargin < 4
        tol = 1e-6;
    end
    if nargin < 5
        max_iter = 10000;
    end

    [n, m] = size(phi_old);  % Get the size of the grid
    sub_dia = deltax_sq * ones(n-3, 1);  % Adjusted length for interior points
    super_dia = sub_dia;
    main_dia = -2 * (deltax_sq + deltay_sq) * ones(n-2, 1);  % Matches the number of interior points in each column

    residuals = zeros(max_iter, 1);  % Array to store residuals at each iteration
    phi_new = phi_old;  % Initialize the new solution

    for iter = 1:max_iter
        % Loop through the columns (excluding boundaries)
        for col = 2:m-1
            % Initialize the right-hand side vector (RHS)
            b = zeros(n-2, 1);

            % Handle the first interior row (i = 2)
            b(1) = source(2, col) * deltax_sq * deltay_sq - phi_old(2, col-1) * deltay_sq - phi_old(2, col+1) * deltay_sq - phi_old(1, col) * deltax_sq;  % Bottom boundary

            % Loop over the interior rows (i = 3 to n-2)
            for row = 3:n-2
                b(row-1) = source(row, col) * deltax_sq * deltay_sq - phi_old(row, col-1) * deltay_sq - phi_old(row, col+1) * deltay_sq;
            end

            % Handle the last interior row (i = n-1)
            b(n-2) = source(n-1, col) * deltax_sq * deltay_sq - phi_old(n-1, col-1) * deltay_sq - phi_old(n-1, col+1) * deltay_sq - phi_old(n, col) * deltax_sq;  % Top boundary

            % Solve the tridiagonal system for this column using TDMA solver
            col_new = tdma_solver(sub_dia, main_dia, super_dia, b);

            % Update the solution for this column
            phi_new(2:n-1, col) = col_new';
        end

        % Calculate the residual as the maximum difference between iterations
        residuals(iter) = max(abs(phi_old(:) - phi_new(:)));

        % Check if convergence is achieved
        if residuals(iter) < tol
            fprintf('Converged in %d iterations with final residual %.6e.\n', iter, residuals(iter));
            residuals = residuals(1:iter);  % Truncate residuals array to actual number of iterations
            num_iters = iter;  % Record the number of iterations
            return;
        end

        % Update old solution for the next iteration
        phi_old = phi_new;
    end

    warning('Maximum iterations reached without convergence.');
    num_iters = max_iter;  % Return max_iter if not converged
end

function solution = tdma_solver(a, b, c, d)
    n = length(d);  % Length of the right-hand side vector

    % Forward sweep
    for i = 2:n
        w = a(i-1) / b(i-1);
        b(i) = b(i) - w * c(i-1);
        d(i) = d(i) - w * d(i-1);
    end

    % Back substitution
    solution = zeros(n, 1);
    solution(n) = d(n) / b(n);

    for i = n-1:-1:1
        solution(i) = (d(i) - c(i) * solution(i+1)) / b(i);
    end
end

