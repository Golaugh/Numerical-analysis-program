% Establish a fixed-point iteration method that satisfies the contraction mapping principle, and calculate the roots of the equation
function q1(N,eps)
% Create the iteration vector and use it as the initial point for iteration
x=zeros(3,1);
x0 = input('Please enter the initial point for the iteration vector (default is 0):');
if x0
    x = x0;
end
% Create the intermediate vector
temp = zeros(3,1);
% Start the iteration
k = 0;
while k < N
    % Use the intermediate vector for communication to ensure synchronization
    temp(1) = cos(x(2)*x(3)+0.5) / 3;
    temp(2) = ((x(1)^2+sin(x(3))+1.06)/81)^0.5 - 0.1;
    temp(3) = (1 - 10*pi/3 - exp(-x(1)*x(2))) / 20;
    % Preemptively check the exit condition for the loop
    min = norm(x-temp,inf);
    % Transfer the intermediate values to the iteration vector
    for j = 1:3
        x(j) = temp(j);
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