clc;
clear all
close all

% Define grid size and spacing
N = 161;
h = 1 / (N-1);

% Initialize the analytical solution field
phi_analytical = zeros(N, N);

% Compute the analytical solution on the grid
for j = 1:N  % y-direction (rows)
    for i = 1:N  % x-direction (columns)
        x = (i-1) * h;
        y = (j-1) * h;
        phi_analytical(j, i) = 500 * exp(-50 * ((1-x)^2 + y^2)) + 100 * x * (1-y);
    end
end

% Generate a meshgrid for plotting, ensuring origin at bottom-left
[X, Y] = meshgrid(0:h:1, 0:h:1);



% Plot the analytical solution
figure;
contourf(X, Y, phi_analytical, 25);
colorbar;
title('Analytical Solution \phi_{analytical}(x, y)');
xlabel('x');
ylabel('y');

