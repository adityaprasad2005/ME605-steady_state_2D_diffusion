function [phi_new, residual_arr] = row_sweep(phi_old, source, deltax_sq, deltay_sq, tol, max_iter)
    if nargin < 4
        tol = 1e-6;
    end
    if nargin < 5
        max_iter = 10000;
    end
    
    residual_arr = zeros(max_iter, 1);  % Array to store residuals for each iteration

    [n, m] = size(phi_old);             % Get the size of the grid
    sub_dia = deltay_sq * ones(m-3, 1);   % Adjusted length for interior points
    super_dia = sub_dia;
    main_dia = -2 * (deltax_sq + deltay_sq) * ones(m-2, 1);  % Matches the number of interior points in each row

    for iter = 1:max_iter
        phi_new = phi_old;

        % Loop through the rows (excluding boundaries)
        for row = 2:n-1
            % Initialize the right-hand side vector (RHS)
            b = zeros(m-2, 1);

            % Handle the first interior column (j = 2)
            b(1) = source(row, 2) * deltax_sq * deltay_sq - phi_old(row-1, 2) * deltax_sq - phi_old(row+1, 2) * deltax_sq - phi_old(row, 1) * deltay_sq;  % Left boundary

            % Loop over the interior columns (j = 3 to m-2)
            for col = 3:m-2
                b(col-1) = source(row, col) * deltax_sq * deltay_sq - phi_old(row-1, col) * deltax_sq - phi_old(row+1, col) * deltax_sq;
            end

            % Handle the last interior column (j = m-1)
            b(m-2) = source(row, m-1) * deltax_sq * deltay_sq - phi_old(row-1, m-1) * deltax_sq - phi_old(row+1, m-1) * deltax_sq - phi_old(row, m) * deltay_sq;  % Right boundary

            % Solve the tridiagonal system for this row using TDMA solver
            row_new = tdma_solver(sub_dia, main_dia, super_dia, b);

            % Update the solution for this row
            phi_new(row, 2:m-1) = row_new';
        end

        % Compute the residual using Residual = A * phi_new - B
        residual = zeros(n, m);
        for i = 2:n-1
            for j = 2:m-1
                % Compute A * phi_new at (i,j)
                A_phi = -deltax_sq * (phi_new(i-1, j) + phi_new(i+1, j)) ...
                        - deltay_sq * (phi_new(i, j-1) + phi_new(i, j+1)) ...
                        + 2 * (deltax_sq + deltay_sq) * phi_new(i, j);

                % Compute B at (i,j)
                B = source(i, j) * deltax_sq * deltay_sq;

                % Residual at (i,j)
                residual(i, j) = A_phi - B;
            end
        end

        % Compute the maximum absolute residual
        residual_arr(iter) = max(abs(residual(:)));

        % Convergence check based on residual
        if residual_arr(iter) < tol
            fprintf('Converged in %d iterations.\n', iter);
            residual_arr = residual_arr(1:iter);  % Trim the residual array to the number of iterations
            phi_new = phi_new;  % Return the updated phi_new
            return;
        end

        % Update old solution for the next iteration
        phi_old = phi_new;
    end

    warning('Maximum iterations reached without convergence.');
    residual_arr = residual_arr(1:max_iter);  % In case of no convergence, return all residuals
end

function solution = tdma_solver(a, b, c, d)
    n = length(d);  % Length of the right-hand side vector

    % Copy b and d to avoid modifying the original vectors
    b_copy = b;
    d_copy = d;

    % Forward sweep
    for i = 2:n
        w = a(i-1) / b_copy(i-1);
        b_copy(i) = b_copy(i) - w * c(i-1);
        d_copy(i) = d_copy(i) - w * d_copy(i-1);
    end

    % Back substitution
    solution = zeros(n, 1);
    solution(n) = d_copy(n) / b_copy(n);

    for i = n-1:-1:1
        solution(i) = (d_copy(i) - c(i) * solution(i+1)) / b_copy(i);
    end
end
