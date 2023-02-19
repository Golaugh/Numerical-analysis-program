%泊松方程列出五点差分格式方程组
clear; close all; clc; 
%默认条件
h1 = 1/10;  j1 = 1/10;  L1 = 1; 
h1_f = input('请输入横向步长:');
j1_f = input('请输入纵向步长:');
if h1_f
    h1 = h1_f;
end
if j1_f
    j1 = j1_f;
end
L_x1 = 1 / h1;  L_t1 = L1 / j1;
%创建右端向量
L_A = L_x1*L_t1;
b = zeros(L_A, 1); 
for i = 1:L_x1
    for j = 1:L_t1
        x_b = i * h1;
        y_b = j * j1;
        b((i-1)*10+j) = (x_b^2 + y_b^2)*exp(x_b*y_b);
    end
end
b = h1^2 * b;
%如需可视化需定义坐标点（目前不需可视化）
%a1 = h1:h1:1;  b1 = j1:j1:L1;
%[x_1, t_1] = meshgrid(a1, b1);
M = zeros(L_A, 1);  
%添加初值条件（边界条件）
for i = 1:L_t1
    M(i,1) = 1;
    M(90+i,1) = exp(i * j1);
    M((i-1)*10 + 1,1) = 1;
    M(i*10,1) = exp(i * h1);
end

%通过推导求得系数矩阵
A = zeros(L_A);  
for i = 1:L_A
    A(i,i) = -4;
    if i-1 > 0
        A(i,i-1) = 1;
    end
    if i+1 < 100
        A(i,i+1) = 1;
    end 
    if i-10 > 0
        A(i,i-10) = 1;
    end 
    if i+10 < 100
        A(i,i+10) = 1;
    end
end
% 显示矩阵A
imagesc(A); 