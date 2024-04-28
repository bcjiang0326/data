function search_combines_betaPtn(n,m,s)
% 2013 0929
% 背景：n 个变量，第 i 个变量水平为 s_i，可以提够0,1, ... , s_i-1 个自由度
% 现在需要 m 个自由度，穷举所有可能的方案。
% 流程介绍：本程序采用递归的方法深度优先检索多叉树
% INPUT:
%           
%           n: 当前待选择的变量个数
%           m: 当前所需自由度的个数
%           s: scalar or n-by-1 vector, 标记各个变量的水平
if n~=round(n) || n < 1
    error('n must be a positive integer!\n');
end
if length(s)==1
    s = s*ones(1,n);
elseif length(s)~=n
    error('s and n do not match!\n');
end
if any(s~=round(s)) || m~=round(m) || min(s)<1 || m<1
    error(' The components of all inputs must be positive integers!\n');
end

b = zeros(1,n);
rec_Ntree(b,n,m,s);

end

function rec_Ntree(b,n,m,s)
%INPUT:   
%   b: 1-by-n vector, 标记各个变量提供的自由度
if sum(s)-n<m
    return;
end
if m == 0
    fprintf('%d ',b);
    fprintf('\n');
    return;
end

for k = s(n)-1:-1:0
    if m >= k
        b(n) = k;
        rec_Ntree(b,n-1,m-k,s(1:n-1));
    end
end
end

%{
% 测试函数1
n = 5; m = 3; s = 3;
search_combines(n,m,s);
% 测试函数2
n = 5; m = 3; s = [3,2,3,2,3];
search_combines(n,m,s);
%}


