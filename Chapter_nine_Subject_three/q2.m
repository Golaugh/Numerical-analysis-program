%主函数(第二问)
%设置基础条件
format long
t0 = 0; tf = 0.004; %t0, tf为区间
Za = [1; 1; 0];     %x初值
h = 0.0001; 
iter = 100;   %迭代次数

%代入梯形方法求解
[t,z] = Trapezoid(@equation, t0, tf , Za, h, iter);

%绘制图形
figure(1)
plot(t,z(:,1),'b',t,z(:,2), 'r',t,z(:,3), 'g--')
title('迭代次数为100')
legend('y1','y2','y3')
figure(2)
plot3(z(:,1),z(:,2),z(:,3));
xlabel('x');ylabel('y');zlabel('z');