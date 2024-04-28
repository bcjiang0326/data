function [Points, NumReps, ID] = checkRep(D,q)
% 20121210
% 检验设计是否存在重复试验点
% 算法: 1.计算出各试验点对应的字典排序法的序号
%       2.将序号按升序排列
%       3.检验是否出现重复的序号并计算重复次数
% INPUT: 
%       D: n-by-s matrix, 被检验的设计
%       q: s-by-1 vector, design各因子水平数
% OUTPUT:
%       Points: t-by-s matrix, 每行代表出现重复的试验点
%       NumReps: t-by-1 vector, Points 中各试验点重复次数
%       ID: t-by-1 vector, Points 中各试验点对应的字典排序的序号。

s = size(q,2);
if s~=1
    error('CheckRep: Input variable q must be a s-by-1 vector!\n');
end
if any( q~=round(q) ) || any( q <= zeros(size(q)) )
    error('Design2Bin:The component of input variable must be positive integer!\n');
end
s = length(q);
if s~= size(D,2)
    error('CheckRep: The sizes of D and q must be matched!\n');
end
n = size(D,1);

% 计算序数
v = D(:,s);
temp = q(s);
for k = s-1:-1:1
    v = v + temp*(D(:,k));
    temp = temp*q(k);
end

v = sort(v);
ID = zeros(0,1);
NumReps = zeros(0,1);
i = 1;
while i < n
    if v(i) == v(i+1)
        ID = cat(1,ID,v(i));
        Num = 2;
        for j = i+2:n
            if v(j) == v(i)
                Num = Num+1;
            else
                break;
            end
        end
        NumReps = cat(1,NumReps,Num);
        i = j;
    else
        i = i+1;
    end
end

if isempty(ID)
    Points = zeros(0,s);
else
    Points = Id2Design(q,ID);
end

end

        


