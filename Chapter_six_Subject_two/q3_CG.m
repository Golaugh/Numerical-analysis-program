%ʹ��CG�����ⷽ����(Matlab�Դ���pcg�����ݲ�����������߽��������ʶ����в�����е�Ԫ�ز�⣬�����ֵ����)
%����ʧ�ܺ��ж�Ϊ�ݶȹ��������Ҫ������ϵ��������ȷ��������������δ�漰��һ������ֵ����������ֵ����ֵ��ϵ
%��˱�ֵ�����ṩ����Ϣ�����޷����뵽�ݶȹ������
function q3_CG(A,b,x0,eps,N)
%Ĭ������
h1 = 1/10;  j1 = 1/10;  L1 = 1; 
L_x1 = 1 / h1;  L_t1 = L1 / j1;
L_A = L_x1*L_t1;
%������ӻ��趨������㣨Ŀǰ������ӻ���
%a1 = h1:h1:1;  b1 = j1:j1:L1;
%[x_1, t_1] = meshgrid(a1, b1);
%�����⺯������
p_x = zeros(L_A,1);
for i = 1:L_x1
    for j = 1:L_t1
        p = (i-1)*10+j;
        p_x(p) = exp(i * h1 * j * j1);
    end
end
%����x��������
x = zeros(L_A, 1);
for i = 1:L_t1
    x(i,1) = 1;
    x(90+i,1) = exp(i * j1);
    x((i-1)*10 + 1,1) = 1;
    x(i*10,1) = exp(i * h1);
end
%����r��������
r = zeros(L_A, N);
%��������������
a = zeros(L_A, 1);
%����beta��������
beta = zeros(L_A, 1);
%����p��������
p = zeros(L_A, 1);
%��p r����������ʼ��
for i = 1:L_A
    p(i) = b(i) - A(i,:)*x0;
    r(i,1) = b(i) - A(i,:)*x0;
end
%���õ�������
count = 0;
%��ʼCG�㷨
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
    %��||x_k+1-x_k||��eps�����㷨ֹͣ�������������ƽ�x_k+1�������������
    min = norm(x-x0,inf);
    if min<eps, break;end
    x0 = x;
end
%�����Ϣ
if count>N
    disp(['��������=  ���㷨����������������',num2str(count)]);
else
    disp(['��������= ',num2str(count)]);
    disp('-------------------------');
    disp(['��⺯�������= ',num2str(norm(x-p_x,inf))]);
end