function aveProjDisc = aveProjDisc_LevelPerm_sym(dim,D,flag)
%20191108 by Bochuan Jiang
%�ο� Jiang and Ai (2019)
%��New construction of uniform designs via average discrepancy theory��
%2.���� D ��level permutation �µ�ƽ�������ԡ�
%   ע�⣺��1��Ҫ�� D ����ͬˮƽ������������
%         ��2��D ���� unbalanced
%INPUT:
%   dim: ͶӰ��ά��
%   D: һ���ǶԳ� OA
%   flag: 'CD','WD','MD'
%OUTPUT:
%   aveProjDisc: ƽ��ƫ��ֵ

n = size(D,2);
s = length(unique(D(:,1)));
if max(size(dim))~=1 || dim~=floor(dim) || dim < 1 || dim > n
    error('Wrong dim!\n')
end
if ~strcmp(flag,'CD') && ~strcmp(flag,'WD') && ~strcmp(flag,'MD')
    error('Wrong flag!\n');
end

%1.����c0��alpha, beta �� gamma
[c0,alpha,beta,gamma] = aveDisc_alpha_beta_gamma(s,flag);
Constant2 = (beta*(gamma+s-1)/s)^dim;
Constant1 = c0^dim-2*alpha^dim+Constant2;
coefvec = zeros(1,dim);
for j = 1:dim
    coefvec(j) = ((gamma-1)/(gamma+s-1))^j * (prod(dim-j+1:dim)/prod(n-j+1:n));
end

%2.���� aveProjDisc
A = GWP(D); %�˺����ڲ�����֤�Գ���
aveProjDisc = Constant1+Constant2*sum(coefvec.*A(1:dim));
end

%{
%���Դ���
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