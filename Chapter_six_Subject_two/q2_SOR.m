% Poisson equation using SOR iteration
function q2_SOR(A,b,x0,eps,w,N)
% Function: Solve the n-order linear system Ax=b using the SOR iteration method

n = length(b); 
x = ones(n,1);
h1 = 1/10;  j1 = 1/10;  L1 = 1; 
L_t1 = L1 / j1;
for i = 1:L_t1
    x(i,1) = 1;
    x(90+i,1) = exp(i * j1);
    x((i-1)*10 + 1,1) = 1;
    x(i*10,1) = exp(i * h1);
end
k = 0;
% Input the coefficient matrix A, right-hand vector b, initial vector x0, tolerance eps, and maximum number of iterations N
% Default conditions
h1 = 1/10;  j1 = 1/10;  L1 = 1; 
L_x1 = 1 / h1;  L_t1 = L1 / j1;
L_A = L_x1 * L_t1;
% Create the solution vector
p_x = zeros(L_A,1);
for i = 1:L_x1
    for j = 1:L_t1
        p = (i-1)*10+j;
        p_x(p) = exp(i * h1 * j * j1);
    end
end
% Perform the iteration while k ≤ N
while k <= N
    % Compute the iterative update for the k-th iteration
    for i = 2:9
        for j = 2:9
            p = (i-1)*10+j;
            x(p) = x0(p) + ( w * (b(p) - A(p,1:p-1) * x(1:p-1) - A(p,p:n) * x0(p:n)) / A(p,p) );
        end
    end
    k = k + 1;
    % If ||x_(k+1) - x_k|| < eps, stop the algorithm and output the approximate solution x_(k+1), otherwise continue iterating
    min = norm(x - x0, inf);
    if min < eps, break; end
    x0 = x;
end
% Output the numerical solution of the system and iteration information
if k > N
    disp(['Number of iterations=  ', num2str(k), ' , algorithm exceeded maximum iterations!']);
else
    disp(['Number of iterations= ', num2str(k)]);
    disp('-------------------------');
    disp(['Error with the solution function= ', num2str(norm(x - p_x, inf))]);
end