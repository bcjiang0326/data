function ptn = MD2_pattern(D,q,dims)
% 2015 1020
% Hickernell and Liu (2002) Uniform designs limit aliasing
% Zhou et al. (2013) mixture discrepancy
% q 水平因子水平集: 0,...,q-1
% 计算投影均匀性
% INPUT:
%       D: N-by-n 设计
%       q: 各因子水平数
%       dims: 1-by-t vector, (k_1,...,k_t), 其中 k_i 为递增正整数，1<=k_i<=n
% OUTPUT:
%       ptn: dims 维投影均匀性

[N,n] = size(D);
if length(q) ~= n || size(q,2)~= 1
    error('q and D do not match each other!\n');
end
if size(dims,1)~=1
    error('dims must be a 1-by-t vector!\n');
end
dims = sort(dims);
if dims(1) < 1 || dims(end)>n || any(dims~=ceil(dims))
    error('k_i must be integers between 1 and s!\n');
end

D = (D+0.5)./(ones(N,1)*(4*q'))-0.125;
ptn = zeros(size(dims));
for kk = 1:length(dims)
    k = dims(kk);
    id = 1:k;
    while id(1)~=-1
        y2 = 0;
        for i = 1:N-1
            for j = i+1:N
                y2 = y2 + prod( 0.875-abs(D(i,id))-abs(D(j,id))-3*abs(D(i,id)-D(j,id))+8*(D(i,id)-D(j,id)).^2 )/N^2;
            end
        end
        y2 = y2*2;

        y2d = 0;
        for i = 1:N
            y2d = y2d + prod( 0.875-2*abs(D(i,id)) )/N^2;
        end

        y1 = 0;
        for i = 1:N
            y1 = y1 + prod( 2/3-abs(D(i,id))-4*D(i,id).^2 )/N;
        end
        y1 = y1*2;

        ptn(kk) = ptn(kk) + (7/12)^k - y1 + y2 + y2d;
        
        id = nchoosek_next(id,n,k);
    end
end
end

%{
% 测试函数
D=[2 3 2; 2 2 2; 3 3 3; 1 1 1; 3 1 1; 3 1 3; 2 2 3;
1 3 3; 2 2 1; 1 1 3; 3 2 2; 2 1 2; 1 3 1; 1 2 2; 3 3 1];
D = D-1; q = [3;3;3]; dims = [1 2 3];
y = MD2_value(D,q);
ptn = MD2_pattern(D,q,dims);
%}