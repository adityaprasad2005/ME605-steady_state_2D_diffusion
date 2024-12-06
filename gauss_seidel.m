function [phi_new, residuals, num_iters] = gauss_seidel(phi_old, source, deltax_sq, deltay_sq, tol, max_iter)
    % Gauss-Seidel method for solving the 2D diffusion equation
    % Inputs:
    %   phi_old    - Initial guess for the solution
    %   source     - Source term (S_phi)
    %   deltax_sq  - Square of grid spacing in x-direction
    %   deltay_sq  - Square of grid spacing in y-direction
    %   tol        - Tolerance for convergence
    %   max_iter   - Maximum number of iterations
    % Outputs:
    %   phi_new    - Final solution after Gauss-Seidel iterations
    %   residuals  - Array of residuals at each iteration
    %   num_iters  - Number of iterations taken to converge

    [n, m] = size(phi_old);  % Get the size of the grid
    phi_new = phi_old;  % Start with the initial guess
    residuals = zeros(max_iter, 1);  % To store residuals at each iteration

    for iter = 1:max_iter
        % Perform Gauss-Seidel update for each interior point
        % Assume all the boundary points to be constant
        for i = 2:n-1
            for j = 2:m-1
                % Update phi at point (i,j) using the finite difference formula
                phi_new(i, j) = (deltay_sq * (phi_old(i+1, j) + phi_old(i-1, j)) + ...
                                deltax_sq * (phi_old(i, j+1) + phi_old(i, j-1)) - ...
                                deltax_sq * deltay_sq * source(i, j)) / (2 * (deltax_sq + deltay_sq));
            end
        end

        % Calculate the residual as the maximum difference between iterations
        residuals(iter) = max(abs(phi_old(:) - phi_new(:)));

        % Check if the solution has converged
        if residuals(iter) < tol
            fprintf('Gauss-Seidel converged in %d iterations with final residual %.6e.\n', iter, residuals(iter));
            residuals = residuals(1:iter);  % Truncate residuals array
            num_iters = iter;  % Record the number of iterations
            return;
        end

        % Update the old solution for the next iteration
        phi_old = phi_new;
    end

    warning('Maximum iterations reached without convergence.');
    num_iters = max_iter;  % Return max_iter if not converged
end

