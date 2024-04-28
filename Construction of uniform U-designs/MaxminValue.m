function y = MaxminValue(L,p)
% 20151015 by xiaopang
% 计算矩阵行距的最小值
% L 为相应的LH, size 为 {N,n}, 元素为 {0,...,N-1}
% p 为采用的范数 
if nargin < 2
    p = 2;
end

epsilon = 1e-12;
N = size(L,1);
if L(1,1)==round(L(1,1)) && L(2,1)==round(L(2,1))
    L = (L+0.5)/N;
end
y = inf;
for i = 1:N-1
    for j = i+1:N
        y1 = norm(L(i,:)-L(j,:),p);
        if y1 < y-epsilon
            y = y1;
        end
    end
end
end