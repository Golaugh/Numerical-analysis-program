n=9;       %网格数
h=1/n;     %步长
F=zeros((n+1)^2,1); %KU=F
K=zeros(size(F));   %K为（n+1）^2阶方阵
U=ones((n+1)^2,1);  %K为（n+1）^2维向量
U0=zeros((n+1)^2,1);%U0用作存边值条件，为初始向量做准备
for i=1:n+1    %对U0中边界点的位置赋值
    U0(i)=1;              %下边界点
    U0(n^2+n+i)=exp(i*h); %右边界点
    U0(1+(i-1)*(n+1))=1;  %左边界点
    U0(i*(n+1))=exp(i*h); %上边界点
end
%% 下面对K赋值
for i=1:(n+1)^2
    K(i,i)=4;
end
for i=1:n-1 %五点中的右点
    for j=2+i*(n+1):(n-1)+i*(n+1)
        K(j,j+1)=-1;
    end
end
for i=1:n-1 %五点中的左点
    for j=3+i*(n+1):n+i*(n+1)
        K(j,j-1)=-1;
    end
end
for i=1:n-2 %五点中的上点
    for j=2+i*(n+1):n+i*(n+1)
        K(j,j+n+1)=-1;
    end
end
for i=2:n-1 %五点中的下点
    for j=2+i*(n+1):n+i*(n+1)
        K(j,j-n-1)=-1;
    end
end
%% 下面对F赋值
%下面赋值非特殊的点
for i=1:n-1
    for j=2+i*(n+1):n+i*(n+1)
        x1 = (floor(j/(n+1)));
        y1 = (mod(j,n+1)-1);
        F(j)=(x1^2+y1^2)*exp(x1*y1*(h^2))*h^4;
    end
end
%下面赋值与边界点相邻的内点(左边和下边为齐次，不用管)
for i=1:n-1 %右边的点
    F(n+i*(n+1))=F(n+i*(n+1))+exp(i*h);  
end
for i=1:n-1 %左边的点
    F(2+i*(n+1))=F(2+i*(n+1))+1;  
end
for i=2:n-1 %下边的点
    F(i+n+1)=F(i+n+1)+1;  
end
for i=1:n-1 %上边的点
    F((n+1)*(n-1)+i+1)=F((n+1)*(n-1)+i+1)+exp(i*h);
end
%下面赋值F中边界点
for i=1:n+1
    F(i)=4;%下边界点
end
for i=1:n-1
    F(1+i*(n+1))=4;%左边界点
    F((i+1)*(n+1))=4*exp(i*h);%右边界点
end
for i=(n+1)*n+1:(n+1)^2 %上边界点
    if i==(n+1)^2
        F(i)=4*exp(n*h);
    else
        F(i)=4*exp((mod(i,n+1)-1)*h);
    end
end

%% 下面对U赋初值
d=sum(sum(U0))/4/n;  %使用边界点值的算术平均值作为U的初值以加快迭代次数
for i=1:(n+1)^2
    U(i)=d;
end
%% 使用高斯赛德迭代求解KU=F
eps=1e-4;  %误差限，可根据需要更改，本次我们使用1e-2和1e-4
N=500;    %最大迭代次数
[u,count] = CG(K,F,U,eps);
%% 下面将解u放到矩阵U1中，并与真解U2作比较
U1=zeros(n+1);%U1用于存放向量u中元素二维化后的数据
U2=zeros(n+1);%U2用于存放真解
for i=1:n+1%u的一维数据放到二维矩阵U1中
    for j=1+(i-1)*(n+1):i*(n+1)
        U1(i,j-(i-1)*(n+1))=u(j);
    end
end
for i=1:n+1 %真解产生的节点函数值放到U2中
    for j=1:n+1
        x2 = (i-1)*h;
        y2 = (j-1)*h;
        U2(i,j)=exp(x2*y2);
    end
end
%% 画图
X=linspace(0,1,n+1);
Y=linspace(0,1,n+1);
mesh(X,Y,U1);
hold on;
title('数值解图像');
%% 创建解函数向量
p_x = zeros(size(F));
for i = 1:n+1
    for j = 1:n+1
        p = (i-1)*10+j;
        p_x(p) = exp(i * h * j * h);
    end
end
%% 输出信息
if count>N
    disp(['迭代次数=  ，算法超出最大迭代次数！',num2str(count)]);
else
    disp(['迭代次数= ',num2str(count)]);
    disp('-------------------------');
    disp('最终向量为= ');
    disp(u);
    disp('-------------------------');
    disp(['与解函数的误差= ',num2str(norm(u-p_x,inf))]);
end