function [cols,m0,uq0,m,uq,D,q] = rand_subOA(D,q0)
%20160519 Xiaopang
%Input: 
%   D: 一个非对称 OA
%   q0: 向量，记录待求子设计各列的水平数
%Output:
%   cols: 子设计列标号
%   m0: 子设计各个水平的因子个数
%   uq0: 子设计出现了哪些水平
%   m: 原设计各个水平的因子个数
%   uq: 原设计出现了哪些水平
%   D：按照水平数由低到高重新排列的原设计
%   q: 重新排列的原设计各个列的水平数

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

uq = unique(q); %各列出现的不同水平数
m = histc(q,[uq-0.5;inf]);
m(end) = []; %等水平因子的个数，m(i) 为水平数为 uq(i) 的因子个数

q0 = sort(q0);
uq0 = unique(q0); %子设计中出现的不同水平数
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

%随机生成一个初始的子设计
cum = cat(1,0,cumsum(m));
cols = zeros(1,0);
for i = 1:length(m0)
    if m0(i)==0
        continue;
    end
    subcols = randperm(m(i),m0(i));
    cols = sort(cat(2,cols,subcols+cum(i)));
end
