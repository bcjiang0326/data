%20160512 by xiaopang
%从特定的若干个集合中抽样，每个集合抽取确定个数的样本，穷举全部的抽样方案
%INPUT:
%   M: 各个元素为各个集合的元素个数
%   M0: 各个元素为各个集合的待抽取的元素个数
%OUTPUT:
%   vec: 每行为一个抽样方案
function vec = enumerate_groups(M,M0)
N = length(M);
id = zeros(1,0);
for i = 1:N
    id = cat(2,id,1:M0(i));
end
vec = zeros(0,length(id));
while id(end)~=-1
    vec = cat(1,vec,id);
    id = choosenextgroup(M,M0,id);
end
end

%{
%测试程序
M = [3;4;5];
M0 = [1;3;4];
vec = enumerate_groups(M,M0);

%}