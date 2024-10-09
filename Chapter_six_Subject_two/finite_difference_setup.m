clear; close all; clc;

% Default parameters for grid resolution and size
h1 = 1/10;  j1 = 1/10;  L1 = 1; 

% Prompt user for step sizes in the horizontal and vertical directions
h1_f = input('Enter horizontal step size: ');
j1_f = input('Enter vertical step size: ');

if h1_f
    h1 = h1_f;
end
if j1_f
    j1 = j1_f;
end

% Calculate the number of grid points in each direction
L_x1 = 1 / h1;  
L_t1 = L1 / j1;
L_A = L_x1 * L_t1;  % Total grid points

% Initialize right-hand side vector based on grid values
b = zeros(L_A, 1); 
for i = 1:L_x1
    for j = 1:L_t1
        x_b = i * h1;
        y_b = j * j1;
        b((i-1)*L_t1 + j) = (x_b^2 + y_b^2) * exp(x_b * y_b);
    end
end
b = h1^2 * b;  % Scale the right-hand side

% Set up boundary conditions in vector M
M = zeros(L_A, 1);  
for i = 1:L_t1
    M(i, 1) = 1;  % Lower boundary
    M((L_x1 - 1) * L_t1 + i, 1) = exp(i * j1);  % Upper boundary
    M((i-1) * L_t1 + 1, 1) = 1;  % Left boundary
    M(i * L_t1, 1) = exp(i * h1);  % Right boundary
end

% Construct coefficient matrix A for the finite difference method
A = zeros(L_A);  
for i = 1:L_A
    A(i,i) = -4;  % Main diagonal
    if i - 1 > 0
        A(i, i - 1) = 1;  % Left neighbor
    end
    if i + 1 <= L_A
        A(i, i + 1) = 1;  % Right neighbor
    end 
    if i - L_t1 > 0
        A(i, i - L_t1) = 1;  % Below neighbor
    end 
    if i + L_t1 <= L_A
        A(i, i + L_t1) = 1;  % Above neighbor
    end
end

% Display the coefficient matrix A visually
imagesc(A);
title('Coefficient Matrix A');
colorbar;
