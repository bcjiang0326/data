function [ID1,D1] = compID(ID0,q)
% 20121211
% 研究补设计的字典排序序号值
% 注：原设计必须无重复
% INPUT:
%       ID0: 原设计的字典排列序号值，必须无重复
%       q: s-by-1 vector, 各个因子水平数
% OUTPUT:
%       ID1: 补设计的字典排列序号值

if size(q,2) ~= 1 || size(ID0,2)~= 1
    error('compID: q and ID0 must be vectors!\n');
end

ID0 = sort(ID0);
N = prod(q);
n = length(ID0);
ID1 = zeros(N-n,1);
k = 1; begin = 0;
for i = 1:n
    for j = begin:ID0(i)-1
        ID1(k) = j;
        k = k+1;
    end
    begin = ID0(i)+1;
end

for j = begin:(N-1)
    ID1(k) = j;
    k = k+1;
end

if nargout > 1
    D1 = Id2Design(q,ID1);
end
end

