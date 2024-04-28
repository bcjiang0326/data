function [cols,m0,uq0,m,uq,D,q] = rand_subOA(D,q0)
%20160519 Xiaopang
%Input: 
%   D: һ���ǶԳ� OA
%   q0: ��������¼��������Ƹ��е�ˮƽ��
%Output:
%   cols: ������б��
%   m0: ����Ƹ���ˮƽ�����Ӹ���
%   uq0: ����Ƴ�������Щˮƽ
%   m: ԭ��Ƹ���ˮƽ�����Ӹ���
%   uq: ԭ��Ƴ�������Щˮƽ
%   D������ˮƽ���ɵ͵����������е�ԭ���
%   q: �������е�ԭ��Ƹ����е�ˮƽ��

n = size(D,2);
q = zeros(n,1);
for i = 1:n
    q(i) = length(unique(D(:,i)));
end
if size(q,2)~=1 || size(q,1)~=n
    error('q must be an n-by-1 vector!\n');
end
if ~issorted(q)
    [q,id] = sort(q);
    D = D(:,id);
end

uq = unique(q); %���г��ֵĲ�ͬˮƽ��
m = histc(q,[uq-0.5;inf]);
m(end) = []; %��ˮƽ���ӵĸ�����m(i) Ϊˮƽ��Ϊ uq(i) �����Ӹ���

q0 = sort(q0);
uq0 = unique(q0); %������г��ֵĲ�ͬˮƽ��
m0 = histc(q0,[uq0-0.5;inf]);
m0(end) = [];

if ~all(ismember(uq0,uq))
    error('The elements of q0 must be the number of certain factor in D !\n');
end

for i = 1:length(uq)
    if ~ismember(uq(i),uq0)
        m0 = [m0(1:i-1);0;m0(i:end)];
    end
end

for i = 1:length(uq)
    if  m0(i) > m(i)
        error('Wrong q0!\n');
    end
end

%�������һ����ʼ�������
cum = cat(1,0,cumsum(m));
cols = zeros(1,0);
for i = 1:length(m0)
    if m0(i)==0
        continue;
    end
    subcols = randperm(m(i),m0(i));
    cols = sort(cat(2,cols,subcols+cum(i)));
end
