% Define the system of equations
function s = equation(t, z)
dy1 = -0.013 * z(1) - 1000 * z(1) * z(2);
dy2 = -2500 * z(2) * z(3);
dy3 = -0.013 * z(1) - 1000 * z(1) * z(2) - 2500 * z(2) * z(3);
s = [dy1 dy2 dy3];   % Return the derivatives as a vector