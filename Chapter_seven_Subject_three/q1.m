%��������ѹ��ӳ��ԭ���Ĳ�����������������㷽�̵ĸ�
function q1(N,eps)
%��������������������Ϊ������ʼ��
x=zeros(3,1);
x0 = input('���������������ʼ�㣨Ĭ��Ϊ0��:');
if x0
    x = x0;
end
%�����м�����
temp = zeros(3,1);
%��ʼ����
k = 0;
while k < N
    %Ϊ��֤ͬ����ʹ���м��������д���
    temp(1) = cos(x(2)*x(3)+0.5) / 3;
    temp(2) = ((x(1)^2+sin(x(3))+1.06)/81)^0.5 - 0.1;
    temp(3) = (1 - 10*pi/3 - exp(-x(1)*x(2))) / 20;
    %��ǰ�ж�ѭ����������
    min = norm(x-temp,inf);
    %���м�ֵ�����������
    for j = 1:3
        x(j) = temp(j);
    end
    %�Ƿ�����ѭ��
    if min<eps, break;end
    k = k + 1;
end
%չʾѭ���������������
if k>N
    disp(['��������=  ���㷨����������������',num2str(k)]);
else
    disp(['��������= ',num2str(k)]);
    disp('-------------------------');
    disp('��÷�����ĸ�Ϊ= ');
    disp(x);
end