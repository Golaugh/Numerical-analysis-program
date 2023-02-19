%泊松方程使用雅可比迭代
function q2_Jacobi(A,b,x0,eps,N)
%功能：用Jacobi迭代法解n阶线性方程组Ax=b

n=length(b); 
x=ones(n,1);
h1 = 1/10;  j1 = 1/10;  L1 = 1; 
L_t1 = L1 / j1;
for i = 1:L_t1
    x(i,1) = 1;
    x(90+i,1) = exp(i * j1);
    x((i-1)*10 + 1,1) = 1;
    x(i*10,1) = exp(i * h1);
end
k=0;
%输入系数矩阵A，右端向量b,以及初始向量x0，精度eps,以及最大迭代次数N
%默认条件
h1 = 1/10;  j1 = 1/10;  L1 = 1; 
L_x1 = 1 / h1;  L_t1 = L1 / j1;
L_A = L_x1*L_t1;
%创建解函数向量
p_x = zeros(L_A,1);
for i = 1:L_x1
    for j = 1:L_t1
        p = (i-1)*10+j;
        p_x(p) = exp(i * h1 * j * j1);
    end
end
%当k≤N时，执行迭代步骤
while k<=N
    %算出第k次迭代递推式
    for i = 2:9
        for j = 2:9
            p = (i-1)*10+j;
            x(p)=(b(p)-A(p,[1:p-1,p+1:n])*x0([1:p-1,p+1:n]))/A(p,p);
        end
    end
    k=k+1;
    %若||x_k+1-x_k||＜eps，则算法停止，输出方程组近似解x_k+1，否则继续迭代
    min = norm(x-x0,inf);
    if min<eps, break;end
    x0 = x;
end
%输出方程组的数值解和迭代信息。
if k>N
    disp(['迭代次数=  ，算法超出最大迭代次数！',num2str(k)]);
else
    disp(['迭代次数= ',num2str(k)]);
    disp('-------------------------');
    disp('最终向量为= ');
    disp(x);
    disp('-------------------------');
    disp(['与解函数的误差= ',num2str(norm(x-p_x,inf))]);
end