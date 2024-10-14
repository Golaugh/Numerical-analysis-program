% Main function
% Set initial conditions
format long
t0 = 0; tf = 0.004; % t0, tf are the interval limits
Za = [1; 1; 0];     % Initial values for x
h = 0.0001; 

% Use the Runge-Kutta method to solve the system
[t, z] = Runge_Kutta(@equation, t0, tf, Za, h);

% Plot the results
figure(1)
plot(t, z(:,1), 'b', t, z(:,2), 'r', t, z(:,3), 'g--')
legend('y1', 'y2', 'y3')
figure(2)
plot3(z(:,1), z(:,2), z(:,3);
xlabel('x'); ylabel('y'); zlabel('z');