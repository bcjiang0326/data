function bptn = beta_pattern_lvl3(D,dims)
% 2013 0927
% 参考 Cheng and Ye (2004)
% Geometric isomorphism and minimum aberration for
% factorial designs with quantitative factors
% INPUT:
%       D: N-by-n 三水平设计，水平集为0,1,2
%       dims: (k_1,...,k_t), 其中 k_i 为递增正整数，1<=k_i<=(3-1)*n
% OUTPUT:
%       ptn: dims 维 beta wordlength pattern
[N,n] = size(D);
dims = sort(dims(:))';
if dims(1) < 1 || dims(end)>2*n || any(dims~=ceil(dims))
    error('k_i must be integers between 1 and s!\n');
end

p1 = [sqrt(1.5),-sqrt(1.5)];
p2 = [1.5*sqrt(2),-3*sqrt(2),0.5*sqrt(2)];

bptn = zeros(size(dims));
for k = 1:length(bptn)
    m = dims(k);
    id = 1:n-1;
    while id(1)~=-1
        order = [id,n+m];
        order(2:end) = order(2:end)-order(1:end-1);
        order = order-1;
        if max(order)>2
            id = nchoosek_next(id,n+m-1,n-1);
            continue;
        end
        Xu = prod(polyval(p1,D(:,order==1)),2)'*prod(polyval(p2,D(:,order==2)),2);
        if abs(Xu)>1e-13
            bptn(k) = bptn(k)+(Xu/N)^2;
        end
        id = nchoosek_next(id,n+m-1,n-1);
    end
end
end




%{
% 测试代码
fcols = 3; q = 3; n = 4; 
N = q^fcols; D0 = zeros(N,n);
D0(:,1:fcols) = combins(fcols,q);
D0(:,4) = mod(2*sum(D0(:,1:fcols)+n,2),3);
D1 = D0; D1(:,4) = mod(2*sum(D1(:,1:fcols),2)+n+1,3);
dims = 4;
bptn0 = beta_pattern_lvl3(D0,dims);
bptn1 = beta_pattern_lvl3(D1,dims);
b0 = 0.5^n*(2+2*(-1)^n);
b1 = 0.5^n*(2 -(-1)^n);
%}