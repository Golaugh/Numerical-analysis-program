%使用牛顿法进行迭代
function q2(N)
%创建迭代向量并将其作为迭代初始点
x=zeros(3,1);
x0 = input('请输入迭代向量初始点（默认为0）:');
if x0
    x = x0;
end
%创建最小误差值
eps=10^-8;
eps0 = input('请输入临界误差值（默认为10^-8）:');
if eps0
    x = eps0;
end
%创建导矩阵
F = zeros(3,3);
%创建中间矩阵
b = zeros(3,1);
%开始迭代
k = 0;
while k < N
    %更新导矩阵
    F(1,1) = 3;
    F(1,2) = x(3)*sin(x(2)*x(3));
    F(1,3) = x(2)*sin(x(2)*x(3));
    F(2,1) = 2*x(1);
    F(2,2) = -162*(x(2)+0.1);
    F(2,3) = cos(x(3));
    F(3,1) = -exp(-x(1)*x(2));
    F(3,2) = -exp(-x(1)*x(2));
    F(3,3) = 20;
    %更新中间矩阵
    b(1) = 3*x(1)-cos(x(2)*x(3))-0.5;
    b(2) = x(1)^2-81*(x(2)+0.1)^2+sin(x(3))+1.06;
    b(3) = exp(-x(1)*x(2))+20*x(3)+10*pi/3-1;
    %更新迭代步伐
    temp = (F^-1) * b;
    %提前判断循环跳出条件
    min = norm(-temp,inf);
    %将中间值传入迭代向量
    for j = 1:3
        x(j) = x(j) - temp(j);
    end
    %是否跳出循环
    if min<eps, break;end
    k = k + 1;
end
%展示循环次数，迭代结果
if k>N
    disp(['迭代次数=  ，算法超出最大迭代次数！',num2str(k)]);
else
    disp(['迭代次数= ',num2str(k)]);
    disp('-------------------------');
    disp('求得方程组的根为= ');
    disp(x);
end