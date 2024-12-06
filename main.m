clc;
clear all
close all

% Define the grid size
N = 21;       % Number of points in x

h_x = 1 / (N-1);  % Grid spacing in x


% Initialize the source term field
S_phi = zeros(N, N);

x = linspace(0,1,N);
y = linspace(0,1,N);


for row_num = 1:1:N         % Row_num refers to row in the Source Matrix. It actually contains the y-coordinates
    for col_num = 1:1:N     % col_num refers to column in the Source Matrix. It actually contains the x-coordinates
        x_dash= x(col_num); % x_dash refers to the x_coordinate in the xy plane
        y_dash = y(row_num);

        S_phi(row_num,col_num) = 50000 * exp(-50*((1-x_dash)^2 + y_dash^2)) * (100*((1-x_dash)^2 + y_dash^2) - 2);
    end
end
% S_phi contains its values which are identically placed as in physical coordinate sense
% S_phi[0] contains the values at y=0 coordinate


% Boundary conditions initialization
boundary_values = zeros(N, N);
for i = 1:N
    x_dash = x(i);
    y_dash = y(i);
    boundary_values(i, 1) = 500 * exp(-50*(1 + y_dash^2));              % Left boundary (x=0)
    boundary_values(i, N) = 100 * (1 - y_dash) + 500 * exp(-50*y_dash^2);    % Right boundary (x=1)
    boundary_values(1, i) = 100 * x_dash + 500 * exp(-50*(1 - x_dash)^2);    % Bottom boundary (y=0)
    boundary_values(N, i) = 500 * exp(-50*((1 - x_dash)^2 + 1));        % Top boundary (y=1)
end




% Initialize coefficient matrix A and right-hand side vector b
A = zeros(N*N, N*N);  % Coefficient matrix
b = zeros(N*N, 1);    % Right-hand side vector


% Construct the coefficient matrix A and vector b
for i = 1:N
    for j = 1:N
        % k refers to node_number
        k = (i-1)*N + j;  % Convert 2D index (i, j) to 1D index k
        
        % COeff Matrix A has k|(i=N, j=N) number of rows and cols
        % Const vector B has k|(i=N, j=N) number of rows 

        if i == 1 || i == N || j == 1 || j == N
            % Boundary points: Apply boundary conditions
            A(k, k) = 1;
            b(k) = boundary_values(i, j);
        else
            % Interior points: Apply the finite difference scheme
            A(k, k) = -4;          % Central point (i, j)
            A(k, k-1) = 1;       % (i-1, j)
            A(k, k+1) = 1;       % (i+1, j)
            A(k, k-N) = 1;       % (i, j-1)
            A(k, k+N) = 1;       % (i, j+1)
            b(k) = (h_x^2 * S_phi(i, j));
        end
    end
end



%%
clc;

tic;  %start time

% Solve the system A*phi = b using Gauss elimination
phi_vector = gaussian_elimination(A, b);

% Reshape the solution vector phi_vector back into a 2D grid
phi = reshape(phi_vector, [N, N]);

elapsed_time = toc; %end time

fprintf('The computational time for solving the PDE is %.6f seconds.\n', elapsed_time);

% Plot the resulting solution
figure();
[X, Y] = meshgrid(0:h_x:1, 0:h_x:1);
contourf(X, Y, phi', 20);
colorbar;
title('Solution \phi(x, y) using Gaussian Elimination for N =', N);
xlabel('x');
ylabel('y');




%%
clc;

%Solve the system of equations A*phi = b using Gauss-seidel iterative method

tic;  %start time
initial_guess = zeros(N, N);
phi_vector = gauss_seidel(boundary_values, S_phi, h_x^2, h_x^2, 1e-6, 100000);

% Reshape the solution vector phi_vector back into a 2D grid
phi = reshape(phi_vector, [N, N]);


elapsed_time = toc;  %end time

fprintf('The computational time for solving the PDE using Gauss-Seidel method is %.6f seconds.\n', elapsed_time);

% Plot the resulting solution
figure();
[X, Y] = meshgrid(0:h_x:1, 0:h_x:1);
contourf(X, Y, phi, 20);
colorbar;
title('Solution \phi(x, y) using Gauss-Seidel iterative method');
xlabel('x');
ylabel('y');

%%
clc;

%Solve the system of equations A*phi = b using row-sweep iterative method

tic;  %start time
[phi_mat, res] = row_sweep(boundary_values, S_phi, h_x^2, h_x^2, 1e-6, 100000);

% Reshape the solution vector phi_vector back into a 2D grid

elapsed_time = toc;  %end time
fprintf('The computational time for solving the PDE using row-sweep method is %.6f seconds.\n', elapsed_time);

% Plot the resulting solution
figure();
[X, Y] = meshgrid(0:h_x:1, 0:h_x:1);
contourf(X, Y, phi_mat, 20);
colorbar;
title('Solution \phi(x, y) using row-sweep iterative method');
xlabel('x');
ylabel('y');

%%
clc;

%Solve the system of equations A*phi = b using column-sweep iterative method

tic; %start time
phi_mat = column_sweep(boundary_values, S_phi, h_x^2, h_x^2, 1e-6, 100000);

% Reshape the solution vector phi_vector back into a 2D grid

elapsed_time = toc;
fprintf('The computational time for solving the PDE using column-sweep method is %.6f seconds.\n', elapsed_time);

% Plot the resulting solution
figure() ;
[X, Y] = meshgrid(0:h_x:1, 0:h_x:1);
contourf(X, Y, phi_mat, 20);
colorbar;
title('Solution \phi(x, y) using column-sweep iterative method');
xlabel('x');
ylabel('y');


%%
clc;

%Solve the system of equations A*phi = b using Alternating Direction method
tic; %start time
phi_mat = ADI(boundary_values, S_phi, h_x^2, h_x^2, 1e-6, 100000);

% Reshape the solution vector phi_vector back into a 2D grid

elapsed_time = toc;
fprintf('The computational time for solving the PDE using ADI method is %.6f seconds.\n', elapsed_time);

% Plot the resulting solution
figure();
[X, Y] = meshgrid(0:h_x:1, 0:h_x:1);
contourf(X, Y, phi_mat, 20);
colorbar;
title('Solution \phi(x, y) using ADI method');
xlabel('x');
ylabel('y');
