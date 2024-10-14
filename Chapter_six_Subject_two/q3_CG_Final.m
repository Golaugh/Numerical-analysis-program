n=9;       % Number of grid points
h=1/n;     % Step size
F=zeros((n+1)^2,1); % KU=F
K=zeros(size(F));   % K is a (n+1)^2 dimensional matrix
U=ones((n+1)^2,1);  % U is a (n+1)^2 dimensional vector
U0=zeros((n+1)^2,1);% U0 is used to store boundary conditions, preparing for the initial vector
for i=1:n+1    % Assign values to the boundary points in U0
    U0(i)=1;              % Lower boundary points
    U0(n^2+n+i)=exp(i*h); % Right boundary points
    U0(1+(i-1)*(n+1))=1;  % Left boundary points
    U0(i*(n+1))=exp(i*h); % Upper boundary points
end

%% Below assigns values to K
for i=1:(n+1)^2
    K(i,i)=4;
end
for i=1:n-1 % Right point of the 5-point stencil
    for j=2+i*(n+1):(n-1)+i*(n+1)
        K(j,j+1)=-1;
    end
end
for i=1:n-1 % Left point of the 5-point stencil
    for j=3+i*(n+1):n+i*(n+1)
        K(j,j-1)=-1;
    end
end
for i=1:n-2 % Upper point of the 5-point stencil
    for j=2+i*(n+1):n+i*(n+1)
        K(j,j+n+1)=-1;
    end
end
for i=2:n-1 % Lower point of the 5-point stencil
    for j=2+i*(n+1):n+i*(n+1)
        K(j,j-n-1)=-1;
    end
end

%% Below assigns values to F
% Assign values to the non-special points
for i=1:n-1
    for j=2+i*(n+1):n+i*(n+1)
        x1 = (floor(j/(n+1)));
        y1 = (mod(j,n+1)-1);
        F(j)=(x1^2+y1^2)*exp(x1*y1*(h^2))*h^4;
    end
end
% Assign values to points adjacent to boundary points (left and bottom are homogeneous, so no need to handle)
for i=1:n-1 % Points on the right
    F(n+i*(n+1))=F(n+i*(n+1))+exp(i*h);  
end
for i=1:n-1 % Points on the left
    F(2+i*(n+1))=F(2+i*(n+1))+1;  
end
for i=2:n-1 % Points at the bottom
    F(i+n+1)=F(i+n+1)+1;  
end
for i=1:n-1 % Points at the top
    F((n+1)*(n-1)+i+1)=F((n+1)*(n-1)+i+1)+exp(i*h);
end
% Assign values to the boundary points in F
for i=1:n+1
    F(i)=4; % Lower boundary points
end
for i=1:n-1
    F(1+i*(n+1))=4; % Left boundary points
    F((i+1)*(n+1))=4*exp(i*h); % Right boundary points
end
for i=(n+1)*n+1:(n+1)^2 % Upper boundary points
    if i==(n+1)^2
        F(i)=4*exp(n*h);
    else
        F(i)=4*exp((mod(i,n+1)-1)*h);
    end
end

%% Below assigns initial values to U
d=sum(sum(U0))/4/n;  % Use the arithmetic mean of boundary points as the initial value for U to speed up the iteration
for i=1:(n+1)^2
    U(i)=d;
end

%% Solve KU=F using Gauss-Seidel iteration
eps=1e-4;  % Error tolerance, can be adjusted as needed; we use 1e-2 and 1e-4 here
N=500;    % Maximum number of iterations
[u,count] = CG(K,F,U,eps);

%% Place the solution u into matrix U1 and compare it with the exact solution U2
U1=zeros(n+1); % U1 stores the 2D version of vector u
U2=zeros(n+1); % U2 stores the exact solution
for i=1:n+1 % Place the 1D data of u into the 2D matrix U1
    for j=1+(i-1)*(n+1):i*(n+1)
        U1(i,j-(i-1)*(n+1))=u(j);
    end
end
for i=1:n+1 % Store the exact solution's node function values into U2
    for j=1:n+1
        x2 = (i-1)*h;
        y2 = (j-1)*h;
        U2(i,j)=exp(x2*y2);
    end
end

%% Plotting
X=linspace(0,1,n+1);
Y=linspace(0,1,n+1);
mesh(X,Y,U1);
hold on;
title('Numerical Solution Plot');

%% Create solution vector
p_x = zeros(size(F));
for i = 1:n+1
    for j = 1:n+1
        p = (i-1)*10+j;
        p_x(p) = exp(i * h * j * h);
    end
end

%% Output information
if count>N
    disp(['Number of iterations=  , algorithm exceeded maximum iterations! ',num2str(count)]);
else
    disp(['Number of iterations= ',num2str(count)]);
    disp('-------------------------');
    disp('Final vector= ');
    disp(u);
    disp('-------------------------');
    disp(['Error with the solution function= ',num2str(norm(u-p_x,inf))]);
end