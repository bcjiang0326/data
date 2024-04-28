function coinDist = CoincidenceDistribution(D,q)
%20190424 by Bochuan Jiang
%���� Coincidence Distribution �ο� Xu(2003) Minimum moment aberration (2)
%ʽ�е�delta_{i,j}(D)
%INPUT:
%   D: N-by-n ��ƣ�Ҫ����ͬˮƽ����������
%   q: n-by-1 vector ��¼ D ���е�ˮƽ����
%OUTPUT:
%   coinDist: N-by-N-by-K vector, K ��ʾ��ˮƽ��ˮƽ����
%       coinDist(i,j,k) ��ʾ��k���е�i�к͵�j����ͬԪ�صĸ���
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





