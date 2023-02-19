%���ɷ���ʹ���ſɱȵ���
function q2_Jacobi(A,b,x0,eps,N)
%���ܣ���Jacobi��������n�����Է�����Ax=b

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
%����ϵ������A���Ҷ�����b,�Լ���ʼ����x0������eps,�Լ�����������N
%Ĭ������
h1 = 1/10;  j1 = 1/10;  L1 = 1; 
L_x1 = 1 / h1;  L_t1 = L1 / j1;
L_A = L_x1*L_t1;
%�����⺯������
p_x = zeros(L_A,1);
for i = 1:L_x1
    for j = 1:L_t1
        p = (i-1)*10+j;
        p_x(p) = exp(i * h1 * j * j1);
    end
end
%��k��Nʱ��ִ�е�������
while k<=N
    %�����k�ε�������ʽ
    for i = 2:9
        for j = 2:9
            p = (i-1)*10+j;
            x(p)=(b(p)-A(p,[1:p-1,p+1:n])*x0([1:p-1,p+1:n]))/A(p,p);
        end
    end
    k=k+1;
    %��||x_k+1-x_k||��eps�����㷨ֹͣ�������������ƽ�x_k+1�������������
    min = norm(x-x0,inf);
    if min<eps, break;end
    x0 = x;
end
%������������ֵ��͵�����Ϣ��
if k>N
    disp(['��������=  ���㷨����������������',num2str(k)]);
else
    disp(['��������= ',num2str(k)]);
    disp('-------------------------');
    disp('��������Ϊ= ');
    disp(x);
    disp('-------------------------');
    disp(['��⺯�������= ',num2str(norm(x-p_x,inf))]);
end