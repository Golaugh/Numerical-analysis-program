function [x, steps] = CG(A, b, x0, eps)
% Conjugate Gradient method for solving A*x = b.
% A - symmetric positive definite matrix.
% b - right-hand side vector.
% x0 - initial guess for the solution.
% eps - tolerance for convergence (optional; default is 1.0e-6).
% x - computed solution.
% steps - number of iterations.

r0 = b - A*x0;  % Initial residual
p0 = r0;        % Initial search direction

if nargin == 3
    eps = 1.0e-6;  % Default tolerance
end

steps = 0;

while true
    if norm(p0) < eps  % Check for convergence
        break;
    end
    
    steps = steps + 1;
    
    a0 = (r0' * r0) / (p0' * A * p0);  % Step size
    x1 = x0 + a0 * p0;                 % Update solution
    r1 = r0 - a0 * A * p0;             % Update residual
    b0 = (r1' * r1) / (r0' * r0);      % Update coefficient
    p1 = r1 + b0 * p0;                 % Update search direction
    
    x0 = x1;  % Prepare for next iteration
    r0 = r1;
    p0 = p1;
end

x = x1;  % Final solution
