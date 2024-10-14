% Use eig to solve the matrix eigenvalues
% Create matrix
A = [10 7 8 7;
    7 5 6 5;
    8 6 10 9;
    7 5 9 10];

B = [2 3 4 5 6;
    4 4 5 6 7;
    0 3 6 7 8;
    0 0 2 8 9;
    0 0 0 1 0];

C = ones(6);
for i = 1:6
    for j = 1:6
        C(i,j) = 1/(i+j-1);
    end
end

% Solve the matrix eigenvalues
a = eig(A);
b = eig(B);
c = eig(C);

% Output information
disp(A); disp("The eigenvalues are:"); disp(a);
disp("----------------------")
disp(B); disp("The eigenvalues are:"); disp(b);
disp("----------------------")
disp(C); disp("The eigenvalues are:"); disp(c);