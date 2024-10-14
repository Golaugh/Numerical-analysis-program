% Solving the system of equations using the CG method (It is not clear how to introduce boundary conditions using Matlab's built-in pcg method, so we decompose the steps from the book and introduce boundary conditions individually)
% After an unsuccessful attempt, it was determined that the main factor influencing the convergence of the conjugate gradient method is the coefficient matrix, and the computation does not yet involve the numerical relationship between one coordinate and neighboring coordinates.
% Therefore, the information from the boundary conditions cannot be incorporated into the conjugate gradient method for now.
function q3_CG(A,b,x0,eps,N)
% Default conditions
h1 = 1/10;  j1 = 1/10;  L1 = 1; 
L_x1 = 1 / h1;  L_t1 = L1 / j1;
L_A = L_x1 * L_t1;
% If visualization is needed, define coordinate points (currently not needed for visualization)
% a1 = h1:h1:1;  b1 = j1:j1:L1;
% [x_1, t_1] = meshgrid(a1, b1);
% Create the solution vector
p_x = zeros(L_A,1);
for i = 1:L_x1
    for j = 1:L_t1
        p = (i-1)*10+j;
        p_x(p) = exp(i * h1 * j * j1);
    end
end
% Create the x-coordinate vector
x = zeros(L_A, 1);
for i = 1:L_t1
    x(i,1) = 1;
    x(90+i,1) = exp(i * j1);
    x((i-1)*10 + 1,1) = 1;
    x(i*10,1) = exp(i * h1);
end
% Create the r iteration vector
r = zeros(L_A, N);
% Create the α iteration vector
a = zeros(L_A, 1);
% Create the beta iteration vector
beta = zeros(L_A, 1);
% Create the p iteration vector
p = zeros(L_A, 1);
% Initialize the p and r iteration vectors
for i = 1:L_A
    p(i) = b(i) - A(i,:) * x0;
    r(i,1) = b(i) - A(i,:) * x0;
end
% Set the number of iterations
count = 0;
% Start the CG algorithm
for i = 1:N

    a(i) = (r(:,i)' * r(:,i)) / (p' * (A * p));

    for i1 = 2:9
        for j1 = 2:9
            pk = (i1-1)*10 + j1;
            x(pk) = x0(pk) + a(i) * p(pk);
        end
    end
    for j = 1:L_A
        r(j,i+1) = r(j,i) - a(i) * (A(j,:) * p);
    end

    beta(i) = (r(:,i+1)' * r(:,i+1)) / (r(:,i)' * r(:,i));
        
    for j = 1:L_A
        p(j) = r(j,i+1) + beta(i) * p(j);
    end
    count = count + 1;
    % If ||x_(k+1) - x_k|| < eps, stop the algorithm and output the approximate solution x_(k+1), otherwise continue iterating
    min = norm(x - x0, inf);
    if min < eps, break; end
    x0 = x;
end
% Output information
if count > N
    disp(['Number of iterations=  ', num2str(count), ' , algorithm exceeded maximum iterations!']);
else
    disp(['Number of iterations= ', num2str(count)]);
    disp('-------------------------');
    disp(['Error with the solution function= ', num2str(norm(x - p_x, inf))]);
end