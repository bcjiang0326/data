function [aveCD,A,BB,VEC] = aveCD_LevelPerm(D,q)
%20181128 by Bochuan Jiang
%�ο� Jiang and Ai (2019) Construction of uniform mixed-level designs
%2.���� D ��level permutation �µ�ƽ��������
%      ��� D �ĸ����Ӳ��ǰ���ˮƽ���������У����������� D �ĸ��У�
%      ʹ�ø���ˮƽ��Ϊ��������
%INPUT:
%   D: һ���ǶԳ� OA
%   q: n-by-1 ��������¼OA���е�ˮƽ��
%OUTPUT:
%   aveCD: ƽ��CDֵ
%   A: generalized wordlength pattern
%   BB: generalized wordtype pattern
%   VEC: ����������������ÿ��Ԫ�صĸ�ˮƽ���ӵĸ���

n = size(D,2);
if nargin < 2 || isempty(q)
    q = zeros(n,1);
    for i = 1:n
        q(i) = length(unique(D(:,i)));
    end
end
if size(q,2)~=1 || size(q,1)~=n
    error('q must be an n-by-1 vector!\n');
end
if ~issorted(q)
    [q,id] = sort(q);
    D = D(:,id);
end


uq = unique(q); %��Ƴ��ֵĲ�ͬˮƽ��ˮƽ������
m = histc(q,[uq-0.5;inf]); %��ˮƽ���Ӹ�������
m(end) = []; %��ˮƽ���ӵĸ�����m(i) Ϊˮƽ��Ϊ uq(i) �����Ӹ���

a0 = zeros(size(uq)); 
a1 = zeros(size(uq)); 
a2 = zeros(size(uq));
for i = 1:length(uq)
    s = uq(i);
    if mod(s,2)==0        
        a0(i) = (26*s^2+1)/(24*s^2);
        a1(i) = (13*s^2+2)/(12*s^2);
        a2(i) = (2*s+2)/(13*s^2+2);
    else
        a0(i) = (13*s^2-1)/(12*s^2);
        a1(i) = (13*s^2-1)/(12*s^2);
        a2(i) = (2*s+2)/(13*s^2-1);
    end
end

BB = genWordTypePtn(D,q);
coef = zeros(1,length(BB));
typeLen = zeros(1,length(BB));
VEC = Id2Design(m+1,(0:length(BB)-1)');

for k = 1:length(BB)
    vec = VEC(k,:);
    coef(k) = prod(a2'.^vec);
    typeLen(k) = sum(vec);
end

A = zeros(1,n);
for j = 1:n
    id = typeLen==j;
    if nargout > 1
        A(j) = sum(BB(id));
    end
end
A(A<1e-12)=0;

aveCD = (13/12)^n-2*prod(a0.^m)+prod(a1.^m)*coef*BB;
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
%}