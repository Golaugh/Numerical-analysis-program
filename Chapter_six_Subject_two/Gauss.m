function [x, k] = Gauss(A, b, x0, ep, N)
% Gauss - Gauss-Seidel iterative method to solve the linear system Ax = b
%
% Inputs:
% A    - Coefficient matrix
% b    - Right-hand side vector
% x0   - Initial guess for the solution (default is the zero vector)
% ep   - Precision (default is 1e-6)
% N    - Maximum number of iterations (default is 5000)
%
% Outputs:
% x    - Solution vector
% k    - Number of iterations taken

n = length(b);

% Set default maximum iterations if N is not provided
if nargin < 5
    N = 5000;
end

% Set default precision if ep is not provided
if nargin < 4
    ep = 1e-6;
end

% Set default initial guess to zero vector if x0 is not provided
if nargin < 3
    x0 = zeros(n, 1);
end

x = zeros(n, 1); % Initialize solution vector
k = 0;           % Initialize iteration counter

while k < N
    for i = 1:n
        if i == 1
            % Update the first element of x
            x(1) = (b(1) - A(1, 2:n) * x0(2:n)) / A(1, 1);
        elseif i == n
            % Update the last element of x
            x(n) = (b(n) - A(n, 1:n-1) * x(1:n-1)) / A(n, n);
        else
            % Update the intermediate elements of x
            x(i) = (b(i) - A(i, 1:i-1) * x(1:i-1) - A(i, i+1:n) * x0(i+1:n)) / A(i, i);
        end
    end
    
    % Check if the solution has converged within the specified tolerance
    if norm(x - x0, inf) < ep
        break;
    end
    
    % Update x0 for the next iteration
    x0 = x;
    k = k + 1;  % Increment iteration counter
end

% Display a warning if the maximum number of iterations is reached
if k == N
    warning('Maximum number of iterations reached');
end

% Display the number of iterations
disp(['k = ', num2str(k)])
end