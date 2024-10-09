function [x, k] = Gauss(A, b, x0, ep, N)
% Gauss-Seidel method for solving the linear system Ax = b.
% A - coefficient matrix.
% b - right-hand side vector.
% x0 - initial guess (optional; default is zero vector).
% ep - tolerance for convergence (optional; default is 1e-6).
% N - maximum number of iterations (optional; default is 5000).
% x - computed solution.
% k - number of iterations.

n = length(b);

if nargin < 5
    N = 5000;  % Default maximum iterations
end
if nargin < 4
    ep = 1e-6;  % Default tolerance
end
if nargin < 3
    x0 = zeros(n, 1);  % Default initial guess
end

x = zeros(n, 1);  % Solution vector
k = 0;

while k < N
    for i = 1:n
        if i == 1
            x(1) = (b(1) - A(1, 2:n) * x0(2:n)) / A(1, 1);
        elseif i == n
            x(n) = (b(n) - A(n, 1:n-1) * x(1:n-1)) / A(n, n);
        else
            x(i) = (b(i) - A(i, 1:i-1) * x(1:i-1) - A(i, i+1:n) * x0(i+1:n)) / A(i, i);
        end
    end
    
    % Check for convergence
    if norm(x - x0, inf) < ep
        break;
    end
    
    x0 = x;  % Update guess for next iteration
    k = k + 1;
end

% Warn if the maximum number of iterations is reached
if k == N
    warning('Maximum number of iterations reached.');
end

disp(['k = ', num2str(k)]);
end
