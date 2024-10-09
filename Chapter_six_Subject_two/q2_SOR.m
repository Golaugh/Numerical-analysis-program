%���ɷ���ʹ��SOR����
function q2_SOR(A,b,x0,eps,w,N)
%���ܣ���SOR��������n�����Է�����Ax=b

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
            x(p)=x0(p) + ( w * (b(p)-A(p,1:p-1)*x(1:p-1)-A(p,p:n)*x0(p:n))/A(p,p) );
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
    disp(['��⺯�������= ',num2str(norm(x-p_x,inf))]);
end