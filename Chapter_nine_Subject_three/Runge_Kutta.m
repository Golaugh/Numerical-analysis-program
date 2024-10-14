% Fourth-order Runge-Kutta method
function [t, z] = Runge_Kutta(fun, t0, tf, Za, h)
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
z(1,:) = Za';             % Align with the direction of variables in the differential equation, change to a row vector

for i = 1:M
    K1 = feval(fun, t(i), z(i,:));                     % K is a row vector
    K2 = feval(fun, t(i) + 1/2 * h, z(i,:) + 1/2 * h * K1);
    K3 = feval(fun, t(i) + 1/2 * h, z(i,:) + 1/2 * h * K2);
    K4 = feval(fun, t(i) + h, z(i,:) + h * K3);   
    z(i + 1, :) = z(i, :) + h / 6 * (K1 + 2 * K2 + 2 * K3 + K4);
end