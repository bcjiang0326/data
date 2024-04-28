function vec = idchangeback(n,k)
%将起点为 1 的下角标 k 转化为2进制, 起点为(0,...0)
%1<=k<=2^n
if k<=0 || k~=round(k) || k > 2^n
    error('k must be a positive integer!\n');
end
vec = zeros(1,n);
if k == 1;
    return;
end
t = 1;
while k-1>=2^t
    t = t+1;
end
vec(t) = 1;
k = k-2^(t-1);
for i = t-1:-1:1
    if k-1>=2^(i-1)
        vec(i) = 1;
        k = k-2^(i-1);
    end
end
vec = fliplr(vec);
if k-1 ~=0
    error('wrong!\n');
end
end
