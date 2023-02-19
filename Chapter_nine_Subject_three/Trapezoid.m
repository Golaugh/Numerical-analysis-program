%梯形方法
function [t,z] = Trapezoid(fun, t0, tf, Za, h, iter)
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
    z(i+1,:) = z(i,:) + h*feval(fun, t(i) , z(i,:));
end

for j = 1:iter
    for i = 1:M
        z(i+1,:) = z(i,:) + h/2 * (feval(fun, t(i) , z(i,:))+feval(fun, t(i+1) , z(i+1,:)));
    end
end