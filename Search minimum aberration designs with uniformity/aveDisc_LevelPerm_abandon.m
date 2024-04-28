function [aveDisc,coinDist] = aveDisc_LevelPerm_abandon(D,q,flag)
%20190416 by Bochuan Jiang
%�ο� Jiang and Ai (2019) Construction of designs with minimum average discrepancy
%2.���� D ��level permutation �µ�ƽ�������ԡ�
%   ע�⣺��1��Ҫ�� D ����ͬˮƽ������������
%         ��2��D ���� unbalanced
%INPUT:
%   D: һ���ǶԳ� OA
%   q: n-by-1 ��������¼OA���е�ˮƽ��
%   flag: 'CD'(default),'WD','MD'
%OUTPUT:
%   aveDisc: ƽ��ƫ��ֵ
%   coinDist: coincidence distributioin 
%       coinDist(i,j,k) ��ʾ��k��group��ˮƽ��ͬ�У��ĵ�i�к�j��Ԫ����ͬ�ĸ���
[N,n] = size(D);
if nargin < 2 || isempty(q)
    q = zeros(n,1);
    for i = 1:n
        q(i) = length(unique(D(:,i)));
    end
end
if ~issorted(q)
    [q,id] = sort(q);
    D = D(:,id);
end
[isgrouped,uq,m] = isgrouped_ByLevNum(q);
if ~isgrouped
    error('q is not grouped by number of levels!\n');
end
if nargin < 3
    flag = 'CD';
end
if ~strcmp(flag,'CD') && ~strcmp(flag,'WD') && ~strcmp(flag,'MD')
    error('Wrong flag!\n');
end


%����coincidence distributioin
coinDist = CoincidenceDistribution(D,q);

%����c0��c1 �� c2
[c0,c1,c2] = aveDisc_c1c2(uq,flag);

%���� alpha
avef1 = zeros(size(uq));
for k = 1:length(uq)
    avef1(k) = sum(kernel_f1(((0:uq(k)-1)'+0.5)/uq(k),flag))/uq(k);
end
alpha = c0^n-2*prod(avef1.^m);

%���� beta
beta = prod((c1./uq./(uq-1)).^m);

Constant = alpha+beta/N*prod(c2.^m);
QuadraticTerm = 0;
for i1 = 1:N-1
    for i2 = i1+1:N
        temp = 1;
        for k = 1:length(uq)
            temp = temp*c2(k)^coinDist(i1,i2,k);
        end
        QuadraticTerm = QuadraticTerm+temp;
    end
end
QuadraticTerm = 2*beta/N^2*QuadraticTerm;
aveDisc = Constant+QuadraticTerm;
end

%{
oa = importdata('OA from web/MA.36.3.7.6.3.finney.txt');
oa = oa(:,5:10);
[N,n] = size(oa);
q = zeros(n,1);
for i = 1:n
    q(i) = length(unique(oa(:,i)));
end
[aveCD,A,BB,VEC] = aveCD_LevelPerm(oa,q);
[aveDisc,coinDist] = aveDisc_LevelPerm_abandon(oa,q,'CD');
%}