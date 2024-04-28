function aveProjDisc = aveProjDisc_LevelPerm_sym(dim,D,flag)
%20191108 by Bochuan Jiang
%参考 Jiang and Ai (2019)
%《New construction of uniform designs via average discrepancy theory》
%2.计算 D 在level permutation 下的平均均匀性。
%   注意：（1）要求 D 中相同水平因子相邻排列
%         （2）D 可以 unbalanced
%INPUT:
%   dim: 投影的维度
%   D: 一个非对称 OA
%   flag: 'CD','WD','MD'
%OUTPUT:
%   aveProjDisc: 平均偏差值

n = size(D,2);
s = length(unique(D(:,1)));
if max(size(dim))~=1 || dim~=floor(dim) || dim < 1 || dim > n
    error('Wrong dim!\n')
end
if ~strcmp(flag,'CD') && ~strcmp(flag,'WD') && ~strcmp(flag,'MD')
    error('Wrong flag!\n');
end

%1.计算c0，alpha, beta 和 gamma
[c0,alpha,beta,gamma] = aveDisc_alpha_beta_gamma(s,flag);
Constant2 = (beta*(gamma+s-1)/s)^dim;
Constant1 = c0^dim-2*alpha^dim+Constant2;
coefvec = zeros(1,dim);
for j = 1:dim
    coefvec(j) = ((gamma-1)/(gamma+s-1))^j * (prod(dim-j+1:dim)/prod(n-j+1:n));
end

%2.计算 aveProjDisc
A = GWP(D); %此函数内部已验证对称性
aveProjDisc = Constant1+Constant2*sum(coefvec.*A(1:dim));
end

%{
%测试代码
s = 3; n = 5; N = 10*s;
q = s*ones(n,1);
D = rand_U_Type_design(q,N);
flag = 'WD';
ProjDisc_ave = zeros(1,n);
ProjDisc = zeros(1,n);
for dim = 1:n
    ProjDisc_ave(dim) = aveProjDisc_LevelPerm_sym(dim,D,flag);
    ProjDisc(dim) = Proj_Disc2(dim, D, flag);
end
clc;
norm(ProjDisc_ave - ProjDisc)
%}