%四阶龙格库塔方法
function [t,z] = Runge_Kutta(fun, t0, tf, Za, h)
%t0, tn为区间
%Za为初值
M = floor((tf-t0)/h) ;      %离散点的个数M+1
if t0 >= tf
    printf('左端点必须小于右端点');
    return;
end
N = length(Za);           %获得变量个数,N
z = zeros(M+1, N);
t =(t0 : h : tf)';
z(1,:) = Za';            %与微分方程中的变量方向统一，变成行向量

for i = 1:M
    K1 =  feval(fun, t(i) , z(i,:));                    %K是行向量
    K2 =  feval(fun, t(i)+1/2*h ,z(i,:)+1/2* h*K1);
    K3 =  feval(fun, t(i)+1/2*h ,z(i,:)+1/2* h*K2);
    K4 =  feval(fun, t(i)+ h ,z(i,:)+ h*K3);   
    z(i+1,:) = z(i,:) + h/6 *(K1 + 2*K2 + 2*K3 + K4);
end