%ʹ��QR�ֽⷨ����������ֵ
function q2_QR(N)
%��������
A = [10 7 8 7;
    7 5 6 5;
    8 6 10 9;
    7 5 9 10];

B = [2 3 4 5 6;
    4 4 5 6 7;
    0 3 6 7 8;
    0 0 2 8 9;
    0 0 0 1 0];

C = ones(6);
for i = 1:6
    for j = 1:6
        C(i,j) = 1/(i+j-1);
    end
end
%���е���
tempA = A;
tempB = B;
tempC = C;
for i = 1:N
    [qa,ra] = qr(tempA);
    tempA = ra * qa;
    
    [qb,rb] = qr(tempB);
    tempB = rb * qb;

    [qc,rc] = qr(tempC);
    tempC = rc * qc;
end
%�������ֵ
eigenA = zeros(length(A),1);
eigenB = zeros(length(B),1);
eigenC = zeros(length(C),1);
for j = 1:length(A)
    eigenA(j) = A(j,j);
end
for j = 1:length(B)
    eigenB(j) = B(j,j);
end
for j = 1:length(C)
    eigenC(j) = C(j,j);
end
%�����Ϣ
disp(A);disp("������ֵΪ:");disp(eigenA);
disp("---------------------------------")
disp(B);disp("������ֵΪ:");disp(eigenB);
disp("---------------------------------")
disp(C);disp("������ֵΪ:");disp(eigenC);