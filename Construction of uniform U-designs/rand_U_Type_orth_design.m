function [U,ID] = rand_U_Type_orth_design(q,n)
% Replication-checking strategy 20130705 
% Obtain a U-type design randomly.
% Input:
%   q: s-vector with positive integer component
%   n: q_1,...,q_s are all divisers of n
% Output:
%   U: U-type design with i-th factor's levels 0,...,q_i-1
%   ID: the index corresponding to each design point

s = length(q);
if min(size(q))~=1
    error('rand_U_Type_orth_design:Input variable must be a vector!\n');
end
r = n./q; % 存储各列的水平重复次数
if any(r~=floor(r))
    error('rand_U_Type_orth_design: n does not match q!\n');
end

% 将 U 按列划分为两个部分，后面分别构造
split = 1;
while split <= s
    if n <= prod(q(1:split))
        break;
    end
    split = split+1;
end


U = zeros(n,s);



%首先构造第一部分
U(:,1)=floor( (0:n-1)'/r(1) );
ID = U(:,1);
for col = 2:split    
    rand_levels = randperm(q(col))-1;
    for i = 1:r(col)
        U((i-1)*q(col)+1:i*q(col),col) = rand_levels;
    end
    if col ~= split
        ID= ID*q(col)+U(:,col);
        [ID,v] = sort(ID);
        U(:,1:col) = U(v,1:col);
    end
end

%构造第二部分
for col = (split+1):s
    U(:,col)=floor( (0:n-1)'/r(col) );
    U(:,col)=U(randperm(n),col);
end
    
if nargout >= 2
    ID = Design2Id(q,U);
end

%{
% 测试程序
levels = 3; s = 5; n = 201;
q = levels*ones(s,1);
 [U,ID] = rand_U_Type_orth_design(q,n);
%}

