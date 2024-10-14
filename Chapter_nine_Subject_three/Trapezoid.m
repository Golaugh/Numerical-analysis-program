% Trapezoidal method
function [t, z] = Trapezoid(fun, t0, tf, Za, h, iter)
% t0, tf are the interval limits
% Za is the initial value
M = floor((tf - t0) / h);      % Number of discrete points is M + 1
if t0 >= tf
    fprintf('The left endpoint must be less than the right endpoint');
    return;
end
N = length(Za);           % Get the number of variables, N
z = zeros(M + 1, N);
t = (t0 : h : tf)';       % Create the time vector
z(1, :) = Za';            % Align with the direction of variables in the differential equation, change to a row vector

for i = 1:M
    z(i + 1, :) = z(i, :) + h * feval(fun, t(i), z(i, :));
end

for j = 1:iter
    for i = 1:M
        z(i + 1, :) = z(i, :) + h / 2 * (feval(fun, t(i), z(i, :)) + feval(fun, t(i + 1), z(i + 1, :)));
    end
end