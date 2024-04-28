function x = randi_SumFixed(isum,imax,m)
%20190412 by Bochuan Jiang
%抽取m个不超过imax 的随机非负整数，使得其和为isum。
x = zeros(1,m);
sum_x = 0;
for i = 1:m
    a = randi(imax+1)-1;
    while sum_x+a > isum || sum_x+a < isum-(m-i)*imax
        a = randi(imax+1)-1;
    end
    x(i) = a;
    sum_x = sum_x+a;
end
end

%{
%测试代码 OA(N,s^n,2)
s = 2; k = 4; N = s^k; 
n_max = (N-1)/(s-1); %饱和设计列数
n = floor(n_max/2);
m = N*(N-1)/2;
M = n*N*(N-s)/(2*s);
imax = (N-s)/s/(s-1); %饱和设计中，任何两行同位置元素相等的个数 
x = randi_SumFixed(M,imax,m);
%}