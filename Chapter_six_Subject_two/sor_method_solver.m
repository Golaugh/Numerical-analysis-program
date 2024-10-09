% Successive Over-Relaxation (SOR) Method
function q2_SOR(A, b, x0, eps, w, N)
% Purpose: Solve Ax = b using the SOR iterative method

n = length(b); 
x = ones(n, 1);
h1 = 1/10; j1 = 1/10; L1 = 1; 
L_t1 = L1 / j1;

% Initialize the solution vector with boundary conditions
for i = 1:L_t1
    x(i, 1) = 1;
    x(90 + i, 1) = exp(i * j1);
    x((i - 1) * 10 + 1, 1) = 1;
    x(i * 10, 1) = exp(i * h1);
end

k = 0;

% Define grid parameters
h1 = 1/10; j1 = 1/10; L1 = 1; 
L_x1 = 1 / h1; L_t1 = L1 / j1;
L_A = L_x1 * L_t1;

% Initialize target solution values
p_x = zeros(L_A, 1);
for i = 1:L_x1
    for j = 1:L_t1
        p = (i - 1) * 10 + j;
        p_x(p) = exp(i * h1 * j * j1);
    end
end

% Iterate until convergence or maximum iterations reached
while k <= N
    % Update the solution vector based on SOR formula
    for i = 2:9
        for j = 2:9
            p = (i - 1) * 10 + j;
            x(p) = x0(p) + (w * (b(p) - A(p, 1:p-1) * x(1:p-1) - A(p, p:n) * x0(p:n)) / A(p, p));
        end
    end
    k = k + 1;

    % Check for convergence
    min = norm(x - x0, inf);
    if min < eps, break; end
    x0 = x;
end

% Output results including iteration count and error norm
if k > N
    disp(['Iteration limit reached: ', num2str(k)]);
else
    disp(['Number of iterations: ', num2str(k)]);
    disp('-------------------------');
    disp(['Maximum error norm: ', num2str(norm(x - p_x, inf))]);
end
