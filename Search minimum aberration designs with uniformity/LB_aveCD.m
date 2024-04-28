function LB = LB_aveCD(N,q)
%20181227 by Bochuan Jiang
%获取全部水平置换下，平均CD值的lower bound
%Input:
%   N: 设计的行数
%   q: n-by-1 vector 设计各列的水平数，要求相同水平因子相邻排列
%Output:
%   LB: Lower Bound

[isgrouped,uq,m] = isgrouped_ByLevNum(q);
if ~isgrouped
    error('LB_aveCD: q is not grouped by number of levels!\n');
end
if any(N./uq~=floor(N./uq))
    error('LB_aveCD: N and q are not match!\n');
end

n = length(q);

alpha1 = 2; %常数项
for k = 1:length(uq)
    if mod(uq(k),2)==0
        alpha1 = alpha1*((26*uq(k)^2+1)/24/uq(k)^2)^m(k);
    else
        alpha1 = alpha1*((13*uq(k)^2-1)/12/uq(k)^2)^m(k);
    end
end
alpha1 = (13/12)^n-alpha1;

c1 = zeros(size(uq));
for k = 1:length(uq)
    if mod(uq(k),2)==0
        c1(k) = (13*uq(k)-2)*(uq(k)-1)/12;
    else
        c1(k) = (13*uq(k)^2-2*uq(k)-3)*(uq(k)-1)/12/uq(k);
    end
end
alpha2 = prod((c1./uq./(uq-1)).^m)/N^2; %二次求和项系数

c2 = zeros(size(uq));
for k = 1:length(uq)
    if mod(uq(k),2)==0
        c2(k) = 15*uq(k)/(13*uq(k)-2);
    else
        c2(k) = (15*uq(k)^2-3)/(13*uq(k)^2-2*uq(k)-3);
    end
end

alpha = alpha1+alpha2*N*prod(c2.^m);

if length(uq)>1
    lb = m.*(N-uq)./(uq*(N-1));
    LB = alpha+alpha2*(N-1)*N*prod(c2.^lb);
else
    a1 = floor(m(k)*(N-uq(k))/(uq(k)*(N-1)));
    r1 = m(k)*N*(N-uq(k))/(2*uq(k))-a1*N*(N-1)/2;
    r2 = N*(N-1)/2-r1;
    LB = alpha+alpha2*2*(r2*c2^a1+r1*c2^(a1+1));
end

%lb = zeros(size(uq));
%for k = 1:length(uq)
    %
%end

end

%{
%测试程序一：Rao-Hamming 方法可达下界
k = 2; s = 5; 
D = RH_OA_prim(k,s);
[N,n] = size(D);
q = s*ones(n,1);
LB = LB_aveCD(N,q);
aveCD = aveCD_LevelPerm(D,q);
LB_ = LB_aveDisc(N,q,'CD');
fprintf('%.8f %.8f %.8f\n',LB,aveCD,LB_);

%测试程序二：MAD.36.6.7.txt，n=s+1,N=s^2, 但s 非素数幂，未达到下界
D = importdata('MAD.36.6/MAD.36.6.7.txt');
[N,n] = size(D);
s = length(unique(D(:,1)));
q = s*ones(n,1);
LB = LB_aveCD(N,q);
aveCD = aveCD_LevelPerm(D,q);
LB_ = LB_aveDisc(N,q,'CD');
fprintf('%.8f %.8f %.8f\n',LB,aveCD,LB_);
[Hamm_pdist,Hamm_Mat] = Hamming_dist(D);

%测试程序三：Rao-Hamming 方法可达下界，混水平
k1 = 2; p1 = 2; m1 = 2; s1 = p1^m1; 
D1 = RH_OA_pow(k1,m1,p1);
n1 = size(D1,2);
k2 = 4; s2 = 2;
D2 = RH_OA_prim(k2,s2);
n2 = size(D2,2);
D = [D1,D2];
[N,n] = size(D);
q = [s1*ones(n1,1);s2*ones(n2,1)];
LB = LB_aveCD(N,q);
aveCD = aveCD_LevelPerm(D,q);
LB_ = LB_aveDisc(N,q,'CD');
fprintf('%.8f %.8f %.8f\n',LB,aveCD,LB_);
%}