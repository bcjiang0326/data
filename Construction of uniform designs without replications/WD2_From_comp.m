function y = WD2_From_comp(q,n,wd)
% 通过 D 的补设计 bar_D 的 WD^2，计算 WD^2(D)
% Input: 
%       q: m-by-1 vector 各元素为各因子的水平数
%       n: D 的试验次数，则 N-n 为 bar_D 的试验次数
%       wd: WD^2(bar_D)
% Output:
%       y: = WD^2(D)
m = length(q);
N = prod(q);

a1 = -(4/3)^m;
a2 = (2*n-N)/n^2;
for k = 1:m
    a2 = a2*(4*q(k)/3+1/(6*q(k)));
end
a3 = (N-n)^2/n^2*(wd+(4/3)^m);
y = a1+a2+a3;

end


