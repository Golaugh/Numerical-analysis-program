function [x, steps] = CG(A, b, x0, eps)
% CG - Conjugate Gradient method for solving Ax = b
%
% Inputs:
% A    - Coefficient matrix
% b    - Right-hand side vector
% x0   - Initial guess for the solution
% eps  - Tolerance for the stopping criterion (optional)
%
% Outputs:
% x     - Solution vector
% steps - Number of iterations taken to converge

r0 = b - A*x0;  % Initial residual
p0 = r0;        % Initial direction is set to the residual

% If only 3 inputs are provided, set the default tolerance
if nargin == 3
    eps = 1.0e-6;
end

steps = 0;

while 1
    % Check the stopping criterion based on the norm of the current direction
    if abs(p0) < eps
        break;
    end
    
    steps = steps + 1;
    
    % Compute step size alpha (r0'*r0 can be reused later, store it)
    a0 = r0' * r0 / (p0' * A * p0);
    
    x1 = x0 + a0 * p0
    
    r1 = r0 - a0 * A * p0;
    
    % Compute beta for the next direction
    b0 = r1' * r1 / (r0' * r0);
    
    % Although r'*r might be reused later, it is not computationally expensive, 
    % so there's no need to store it as a new variable.
    
    % Update direction
    p1 = r1 + b0 * p0;
    
    % Since only the previous and current vectors are needed, we can overwrite them
    x0 = x1;
    r0 = r1;
    p0 = p1;
end

% Return the final solution
x = x0;
end
