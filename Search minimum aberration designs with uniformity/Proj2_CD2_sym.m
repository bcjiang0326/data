function y = Proj2_CD2_sym(D)
% 20200211 参考《Uniform Projection designs》 Theorem 2
% INPUT:
%       D: 设计   
% OUTPUT:
%       y: CD准则下的二维投影均匀性

[N,n] = size(D);
s = length(unique(D(:,1)));
for j = 2:n
    if s ~= length(unique(D(:,j)))
        error('D must be symmetrical!\n');
    end
end
pd = squareform(pdist(D,'minkowski',1));
C = (4*(5*n-2)*s^4+30*(3*n-5)*s^2+15*n+33)/(720*(n-1)*s^4);
C = C + (1+(-1)^s)/(64*s^4);

g1 = 0;
for i = 1:N-1
    for j = i+1:N
        g1 = g1+pd(i,j)^2;
    end
end
g1 = g1*2;

g2 = 0;
for i = 1:N
    g2 = g2+(sum(pd(i,:)))^2;
end
g2 = g2*2/N;

g = g1-g2;

y = g/(4*n*(n-1)*N^2*s^2)+C;
end
%{
%测试代码
s = 4; n = 5; N = 3*s;
q = s*ones(n,1);
D = rand_U_Type_design(q,N);
flag = 'CD';
Proj_Disc = Proj_Disc2(2, D, flag);
Proj_CD2 = Proj2_CD2_sym(D);
fprintf('%.8f  %.8f\n',Proj_Disc,Proj_CD2); 
%}

