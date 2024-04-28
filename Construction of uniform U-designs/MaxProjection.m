function y = MaxProjection(L)
% 20151019 by xiaopang
% 参考 Maximum projection designs for computer experiments 中定理 1
% L 为 Latin hypercube，元素为 0,...,N-1
[N,n] = size (L);
L = (L+0.5)/N;
y = 0;
for i = 1:N-1
    for j = i+1:N
        y = y+1/(prod(L(i,:)-L(j,:))^2);
    end
end
m = nchoosek(N,2);
y = (y/m)^(1/n);