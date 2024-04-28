function B = dist_distr_multi(D,q)
%20150616 by Xiaopang
%计算 D 的 multiple distance distribution
%INPUT:
%   D: 一个非对称 OA
%   q: n-by-1 向量，记录OA各列的水平数
%OUTPUT:
%   B: 计算 D 的 multiple distance distribution
[N,n] = size(D);
if size(q,2)~=1 || size(q,1)~=n
    error('q must be an n-by-1 vector!\n');
end
B = zeros(2^n,1);
for i = 1:N-1
    for j = i+1:N
        vec = (D(i,:)-D(j,:))~=0;
        k = polyval(vec,2)+1;
        B(k) = B(k)+1;
    end
end
B = B*2;
B(1) = B(1)+N;
B = B/N;
end




