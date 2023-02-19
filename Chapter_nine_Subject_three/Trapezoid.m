%���η���
function [t,z] = Trapezoid(fun, t0, tf, Za, h, iter)
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
    z(i+1,:) = z(i,:) + h*feval(fun, t(i) , z(i,:));
end

for j = 1:iter
    for i = 1:M
        z(i+1,:) = z(i,:) + h/2 * (feval(fun, t(i) , z(i,:))+feval(fun, t(i+1) , z(i+1,:)));
    end
end