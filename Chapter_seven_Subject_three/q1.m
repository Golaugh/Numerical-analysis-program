%建立满足压缩映射原理的不动点迭代法，并计算方程的根
function q1(N,eps)
%创建迭代向量并将其作为迭代初始点
x=zeros(3,1);
x0 = input('请输入迭代向量初始点（默认为0）:');
if x0
    x = x0;
end
%创建中间向量
temp = zeros(3,1);
%开始迭代
k = 0;
while k < N
    %为保证同步性使用中间向量进行传递
    temp(1) = cos(x(2)*x(3)+0.5) / 3;
    temp(2) = ((x(1)^2+sin(x(3))+1.06)/81)^0.5 - 0.1;
    temp(3) = (1 - 10*pi/3 - exp(-x(1)*x(2))) / 20;
    %提前判断循环跳出条件
    min = norm(x-temp,inf);
    %将中间值传入迭代向量
    for j = 1:3
        x(j) = temp(j);
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