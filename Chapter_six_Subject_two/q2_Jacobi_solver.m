function q2_Jacobi(A, b, x0, eps, N)
% Function: Solves the linear system Ax = b using the Jacobi iterative method.
% Inputs:
% A  - Coefficient matrix.
% b  - Right-hand side vector.
% x0 - Initial guess for the solution.
% eps - Tolerance for convergence (optional; default is 1e-6).
% N  - Maximum number of iterations (optional; default is 500).

n = length(b);  % Get the size of the system
x = ones(n, 1);  % Initial solution vector

% Boundary conditions setup
h1 = 1/10;  j1 = 1/10;  L1 = 1;
L_t1 = L1 / j1;
for i = 1:L_t1
    x(i, 1) = 1;
    x(90 + i, 1) = exp(i * j1);
    x((i - 1) * 10 + 1, 1) = 1;
    x(i * 10, 1) = exp(i * h1);
end

% Default tolerance and max iterations if not provided
if nargin < 5
    N = 500;  % Default max iterations
end
if nargin < 4
    eps = 1e-6;  % Default tolerance
end
if nargin < 3
    x0 = zeros(n, 1);  % Default initial guess
end

k = 0;  % Initialize iteration counter

% Precompute exact solution (if necessary)
L_x1 = 1 / h1;  
L_A = L_x1 * L_t1;
p_x = zeros(L_A, 1);
for i = 1:L_x1
    for j = 1:L_t1
        p = (i - 1) * 10 + j;
        p_x(p) = exp(i * h1 * j * j1);
    end
end

% Jacobi iteration loop
while k <= N
    % Compute next iteration
    for i = 2:9
        for j = 2:9
            p = (i - 1) * 10 + j;
            x(p) = (b(p) - A(p, [1:p-1, p+1:n]) * x0([1:p-1, p+1:n])) / A(p, p);
        end
    end
    k = k + 1;
    
    % Check for convergence
    if norm(x - x0, inf) < eps
        break;
    end
    x0 = x;  % Update solution for next iteration
end

% Display results
if k > N
    disp(['Number of iterations = ', num2str(k), ' (Reached max iteration limit).']);
else
    disp(['Number of iterations = ', num2str(k)]);
    disp('Final solution:');
    disp(x);
    disp('-------------------------');
    disp(['Error with exact solution = ', num2str(norm(x - p_x, inf))]);
end
