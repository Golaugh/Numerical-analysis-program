n=9;       %������
h=1/n;     %����
F=zeros((n+1)^2,1); %KU=F
K=zeros(size(F));   %KΪ��n+1��^2�׷���
U=ones((n+1)^2,1);  %KΪ��n+1��^2ά����
U0=zeros((n+1)^2,1);%U0�������ֵ������Ϊ��ʼ������׼��
for i=1:n+1    %��U0�б߽���λ�ø�ֵ
    U0(i)=1;              %�±߽��
    U0(n^2+n+i)=exp(i*h); %�ұ߽��
    U0(1+(i-1)*(n+1))=1;  %��߽��
    U0(i*(n+1))=exp(i*h); %�ϱ߽��
end
%% �����K��ֵ
for i=1:(n+1)^2
    K(i,i)=4;
end
for i=1:n-1 %����е��ҵ�
    for j=2+i*(n+1):(n-1)+i*(n+1)
        K(j,j+1)=-1;
    end
end
for i=1:n-1 %����е����
    for j=3+i*(n+1):n+i*(n+1)
        K(j,j-1)=-1;
    end
end
for i=1:n-2 %����е��ϵ�
    for j=2+i*(n+1):n+i*(n+1)
        K(j,j+n+1)=-1;
    end
end
for i=2:n-1 %����е��µ�
    for j=2+i*(n+1):n+i*(n+1)
        K(j,j-n-1)=-1;
    end
end
%% �����F��ֵ
%���渳ֵ������ĵ�
for i=1:n-1
    for j=2+i*(n+1):n+i*(n+1)
        x1 = (floor(j/(n+1)));
        y1 = (mod(j,n+1)-1);
        F(j)=(x1^2+y1^2)*exp(x1*y1*(h^2))*h^4;
    end
end
%���渳ֵ��߽�����ڵ��ڵ�(��ߺ��±�Ϊ��Σ����ù�)
for i=1:n-1 %�ұߵĵ�
    F(n+i*(n+1))=F(n+i*(n+1))+exp(i*h);  
end
for i=1:n-1 %��ߵĵ�
    F(2+i*(n+1))=F(2+i*(n+1))+1;  
end
for i=2:n-1 %�±ߵĵ�
    F(i+n+1)=F(i+n+1)+1;  
end
for i=1:n-1 %�ϱߵĵ�
    F((n+1)*(n-1)+i+1)=F((n+1)*(n-1)+i+1)+exp(i*h);
end
%���渳ֵF�б߽��
for i=1:n+1
    F(i)=4;%�±߽��
end
for i=1:n-1
    F(1+i*(n+1))=4;%��߽��
    F((i+1)*(n+1))=4*exp(i*h);%�ұ߽��
end
for i=(n+1)*n+1:(n+1)^2 %�ϱ߽��
    if i==(n+1)^2
        F(i)=4*exp(n*h);
    else
        F(i)=4*exp((mod(i,n+1)-1)*h);
    end
end

%% �����U����ֵ
d=sum(sum(U0))/4/n;  %ʹ�ñ߽��ֵ������ƽ��ֵ��ΪU�ĳ�ֵ�Լӿ��������
for i=1:(n+1)^2
    U(i)=d;
end
%% ʹ�ø�˹���µ������KU=F
eps=1e-4;  %����ޣ��ɸ�����Ҫ���ģ���������ʹ��1e-2��1e-4
N=500;    %����������
[u,count] = CG(K,F,U,eps);
%% ���潫��u�ŵ�����U1�У��������U2���Ƚ�
U1=zeros(n+1);%U1���ڴ������u��Ԫ�ض�ά���������
U2=zeros(n+1);%U2���ڴ�����
for i=1:n+1%u��һά���ݷŵ���ά����U1��
    for j=1+(i-1)*(n+1):i*(n+1)
        U1(i,j-(i-1)*(n+1))=u(j);
    end
end
for i=1:n+1 %�������Ľڵ㺯��ֵ�ŵ�U2��
    for j=1:n+1
        x2 = (i-1)*h;
        y2 = (j-1)*h;
        U2(i,j)=exp(x2*y2);
    end
end
%% ��ͼ
X=linspace(0,1,n+1);
Y=linspace(0,1,n+1);
mesh(X,Y,U1);
hold on;
title('��ֵ��ͼ��');
%% �����⺯������
p_x = zeros(size(F));
for i = 1:n+1
    for j = 1:n+1
        p = (i-1)*10+j;
        p_x(p) = exp(i * h * j * h);
    end
end
%% �����Ϣ
if count>N
    disp(['��������=  ���㷨����������������',num2str(count)]);
else
    disp(['��������= ',num2str(count)]);
    disp('-------------------------');
    disp('��������Ϊ= ');
    disp(u);
    disp('-------------------------');
    disp(['��⺯�������= ',num2str(norm(u-p_x,inf))]);
end