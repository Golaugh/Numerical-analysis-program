%使用CG方法解方程组(Matlab自带的pcg方法暂不清楚如何引入边界条件，故对书中步骤进行单元素拆解，引入边值条件)
%尝试失败后判断为梯度共轭法迭代主要因素由系数矩阵所确定，且运算中尚未涉及到一个坐标值与相邻坐标值的数值关系
%因此边值条件提供的信息尚且无法加入到梯度共轭方法中
function q3_CG(A,b,x0,eps,N)
%默认条件
h1 = 1/10;  j1 = 1/10;  L1 = 1; 
L_x1 = 1 / h1;  L_t1 = L1 / j1;
L_A = L_x1*L_t1;
%如需可视化需定义坐标点（目前不需可视化）
%a1 = h1:h1:1;  b1 = j1:j1:L1;
%[x_1, t_1] = meshgrid(a1, b1);
%创建解函数向量
p_x = zeros(L_A,1);
for i = 1:L_x1
    for j = 1:L_t1
        p = (i-1)*10+j;
        p_x(p) = exp(i * h1 * j * j1);
    end
end
%创建x坐标向量
x = zeros(L_A, 1);
for i = 1:L_t1
    x(i,1) = 1;
    x(90+i,1) = exp(i * j1);
    x((i-1)*10 + 1,1) = 1;
    x(i*10,1) = exp(i * h1);
end
%创建r迭代向量
r = zeros(L_A, N);
%创建α迭代向量
a = zeros(L_A, 1);
%创建beta迭代向量
beta = zeros(L_A, 1);
%创建p迭代向量
p = zeros(L_A, 1);
%将p r迭代向量初始化
for i = 1:L_A
    p(i) = b(i) - A(i,:)*x0;
    r(i,1) = b(i) - A(i,:)*x0;
end
%设置迭代次数
count = 0;
%开始CG算法
for i = 1:N

    a(i) = (r(:,i)' * r(:,i)) / (p' * (A * p));

    for i1 = 2:9
        for j1 = 2:9
            pk = (i1-1)*10+j1;
            x(pk) = x0(pk) + a(i) * p(pk);
        end
    end
    for j = 1:L_A
        r(j,i+1) = r(j,i) - a(i) *  (A(j,:) * p);
    end

    beta(i) = (r(:,i+1)'* r(:,i+1)) / (r(:,i)' * r(:,i));
        
    for j = 1:L_A
        p(j) = r(j,i+1) + beta(i) * p(j);
    end
    count = count + 1;
    %若||x_k+1-x_k||＜eps，则算法停止，输出方程组近似解x_k+1，否则继续迭代
    min = norm(x-x0,inf);
    if min<eps, break;end
    x0 = x;
end
%输出信息
if count>N
    disp(['迭代次数=  ，算法超出最大迭代次数！',num2str(count)]);
else
    disp(['迭代次数= ',num2str(count)]);
    disp('-------------------------');
    disp(['与解函数的误差= ',num2str(norm(x-p_x,inf))]);
end