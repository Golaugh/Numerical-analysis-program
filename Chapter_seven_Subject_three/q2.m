% Use Newton's method for iteration
function q2(N)
% Create the iteration vector and use it as the initial point for iteration
x=zeros(3,1);
x0 = input('Please enter the initial point for the iteration vector (default is 0):');
if x0
    x = x0;
end
% Create the minimum error value
eps=10^-8;
eps0 = input('Please enter the critical error value (default is 10^-8):');
if eps0
    eps = eps0;
end
% Create the Jacobian matrix
F = zeros(3,3);
% Create the intermediate vector
b = zeros(3,1);
% Start the iteration
k = 0;
while k < N
    % Update the Jacobian matrix
    F(1,1) = 3;
    F(1,2) = x(3)*sin(x(2)*x(3));
    F(1,3) = x(2)*sin(x(2)*x(3));
    F(2,1) = 2*x(1);
    F(2,2) = -162*(x(2)+0.1);
    F(2,3) = cos(x(3));
    F(3,1) = -exp(-x(1)*x(2));
    F(3,2) = -exp(-x(1)*x(2));
    F(3,3) = 20;
    
    % Update the intermediate vector
    b(1) = 3*x(1)-cos(x(2)*x(3))-0.5;
    b(2) = x(1)^2-81*(x(2)+0.1)^2+sin(x(3))+1.06;
    b(3) = exp(-x(1)*x(2))+20*x(3)+10*pi/3-1;
    
    % Update the iteration step
    temp = (F^-1) * b;
    
    % Preemptively check the exit condition for the loop
    min = norm(-temp,inf);
    % Transfer the intermediate values to the iteration vector
    for j = 1:3
        x(j) = x(j) - temp(j);
    end
    
    % Check if the loop should exit
    if min<eps, break;end
    k = k + 1;
end

% Display the number of iterations and the results
if k>N
    disp(['Number of iterations=  , algorithm exceeded maximum iterations! ',num2str(k)]);
else
    disp(['Number of iterations= ',num2str(k)]);
    disp('-------------------------');
    disp('The roots of the equation system are= ');
    disp(x);
end