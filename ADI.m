function [phi_new, residuals, num_iters] = ADI(phi_old, source, deltax_sq, deltay_sq, tol, max_iter)
    % ADI_method implements the ADI method using row and column sweeps.
    % Inputs:
    %   phi_old    - Initial guess for the solution.
    %   source     - Source term (S_phi).
    %   deltax_sq  - Square of grid spacing in the x-direction.
    %   deltay_sq  - Square of grid spacing in the y-direction.
    %   tol        - Tolerance for convergence.
    %   max_iter   - Maximum number of iterations.
    % Outputs:
    %   phi_new    - Final solution after ADI iterations.
    %   residuals  - Array of residuals at each iteration.
    %   num_iters  - Number of iterations taken to converge.

    % Initialize
    [n, m] = size(phi_old);
    phi_new = phi_old;  % Start with the initial guess
    residuals = zeros(max_iter, 1);  % To store residuals at each iteration

    for iter = 1:max_iter
        % Perform row sweep
        phi_mid = row_sweep(phi_new, source, deltax_sq, deltay_sq, tol, 1);

        % Perform column sweep
        phi_new = column_sweep(phi_mid, source, deltax_sq, deltay_sq, tol, 1);

        % Calculate the residual as the maximum difference between iterations
        residuals(iter) = max(abs(phi_old(:) - phi_new(:)));

        % Check for convergence
        if residuals(iter) < tol
            %fprintf('ADI converged in %d iterations with final residual %.6e.\n', iter, residuals(iter));
            residuals = residuals(1:iter);  % Truncate residuals array
            num_iters = iter;  % Record the number of iterations
            return ;
        end

        %fprintf('ADI method converged in %d iterations with final residual')

        % Update the old solution for the next iteration
        phi_old = phi_new;
    end

    % warning('Maximum iterations reached without convergence.');
    num_iters = 2*max_iter;  % Return max_iter if not converged
end

