function [h,b,H,B] = WD_GWPoa_lhd(N,n,s)
%20150603 by Xiaopang
%aveWD = h*A(D)+b
%aveWD_j = H(i,:)*A(D)+B(i)
%Let PWD = (aveWD_0,,...,aveWD_n)', A(D) = (A_0(D),A_1(D),...,A_n(D))'
%Then PWD = H*A(D)+B
%Input:
%       N: 行数
%       n: 列数
%       s: OA 水平数

b = -(4/3)^n+1.5^n/N-(1.5-1/(3*s)-1/(3*N)+1/(6*s^2)+1/(6*N*s))^n/N;
a1 = 4/3-1/(3*N*s)-1/(6*N^2)+1/(6*N*s^2)+1/(6*N^2*s);
a2 = (N^2*s-N^2-2*N*s+N+s)/(8*N^2*s^2-2*N*s-s^2+N+s);
h = zeros(1,n+1);
for i = 1:n+1
    h(i) = a1^n*a2^(i-1);
end

if nargout > 2
    B = zeros(n+1,1);
    for j = 1:n+1
        %B(j) = -nchoosek(n,j-1)*(1/3)^(j-1);
        B(j) = nchoosek(n,j-1)*...
            (-(1/3)^(j-1)+0.5^(j-1)/N-(0.5-1/(3*s)-1/(3*N)+1/(6*s^2)+1/(6*N*s))^(j-1)/N);
    end
    a1hat = a1-1;
    a2hat = (N^2*s-N^2-2*N*s+N+s)/(2*N^2*s^2-2*N*s-s^2+N+s);
    H = zeros(n,n+1);
    for j = 1:n+1
        for i = 1:j
            H(j,i) = a1hat^(j-1)*a2hat^(i-1)*nchoosek(n-i+1,j-i);
        end
    end
end

end
%{ 
%测试代码
N = 12; s = 3; n = 4;
[h,b,H,B] = WD_GWPoa_lhd(N,n,s);
z1 = norm(sum(H)-h);
z2 = norm(sum(B)-b);
%}