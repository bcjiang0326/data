function coinDist = CoincidenceDistribution(D,q)
%20190424 by Bochuan Jiang
%计算 Coincidence Distribution 参考 Xu(2003) Minimum moment aberration (2)
%式中的delta_{i,j}(D)
%INPUT:
%   D: N-by-n 设计，要求相同水平数的列相邻
%   q: n-by-1 vector 记录 D 各列的水平数。
%OUTPUT:
%   coinDist: N-by-N-by-K vector, K 表示混水平的水平组数
%       coinDist(i,j,k) 表示第k组中第i行和第j行相同元素的个数
[N,n] = size(D);
if nargin < 2 || isempty(q)
    q = zeros(n,1);
    for i = 1:n
        q(i) = length(unique(D(:,i)));
    end
end
[isgrouped,uq,m] = isgrouped_ByLevNum(q);
if ~isgrouped
    error('q is not grouped by number of levels!\n');
end

coinDist = zeros(N,N,length(uq));
ed = 0;
for k = 1:length(uq)
    ed = ed+m(k);
    bg = ed-m(k)+1;
    for i = 1:N
        coinDist(i,i,k) = m(k);
        for j = i+1:N
            coinDist(i,j,k) = sum(D(i,bg:ed)==D(j,bg:ed));
            coinDist(j,i,k) = coinDist(i,j,k);
        end
    end
end
end





