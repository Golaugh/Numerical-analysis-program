%������(�ڶ���)
%���û�������
format long
t0 = 0; tf = 0.004; %t0, tfΪ����
Za = [1; 1; 0];     %x��ֵ
h = 0.0001; 
iter = 100;   %��������

%�������η������
[t,z] = Trapezoid(@equation, t0, tf , Za, h, iter);

%����ͼ��
figure(1)
plot(t,z(:,1),'b',t,z(:,2), 'r',t,z(:,3), 'g--')
title('��������Ϊ100')
legend('y1','y2','y3')
figure(2)
plot3(z(:,1),z(:,2),z(:,3));
xlabel('x');ylabel('y');zlabel('z');