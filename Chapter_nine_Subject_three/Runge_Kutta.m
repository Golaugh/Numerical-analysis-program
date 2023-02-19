%�Ľ������������
function [t,z] = Runge_Kutta(fun, t0, tf, Za, h)
%t0, tnΪ����
%ZaΪ��ֵ
M = floor((tf-t0)/h) ;      %��ɢ��ĸ���M+1
if t0 >= tf
    printf('��˵����С���Ҷ˵�');
    return;
end
N = length(Za);           %��ñ�������,N
z = zeros(M+1, N);
t =(t0 : h : tf)';
z(1,:) = Za';            %��΢�ַ����еı�������ͳһ�����������

for i = 1:M
    K1 =  feval(fun, t(i) , z(i,:));                    %K��������
    K2 =  feval(fun, t(i)+1/2*h ,z(i,:)+1/2* h*K1);
    K3 =  feval(fun, t(i)+1/2*h ,z(i,:)+1/2* h*K2);
    K4 =  feval(fun, t(i)+ h ,z(i,:)+ h*K3);   
    z(i+1,:) = z(i,:) + h/6 *(K1 + 2*K2 + 2*K3 + K4);
end