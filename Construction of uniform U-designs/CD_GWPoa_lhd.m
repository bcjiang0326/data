function [h,b,H,B] = CD_GWPoa_lhd(N,n,s)
%20150603 by Xiaopang
%aveCD = h*A(D)+b
%aveCD_j = H(i,:)*A(D)+B(i)
%Let PCD = (aveCD_0,,...,aveCD_n)', A(D) = (A_0(D),A_1(D),...,A_n(D))'
%Then PCD = H*A(D)+B
%Input:
%       N: 行数
%       n: 列数
%       s: OA 水平数

if mod(N,2)==0
    b = (13/12)^n-2*(13/12+1/(24*N^2))^n+1.25^n/N-(1.25-1/(6*s)-1/(6*N))^n/N;
    a1 = 13/12-1/(6*N*s);
    a2 = (2*N-2)/(13*N*s-2);
else
    b = (13/12)^n-2*(13/12-1/(12*N^2))^n+(1.25-1/(4*N^2))^n/N-(1.25-1/(6*s)-1/(6*N)-1/(4*N^2))^n/N;
    a1 = 13/12-1/(6*N*s)-1/(4*N^2);
    a2 = (2*N-2)/(13*N*s-2-3*s/N);
end
h = zeros(1,n+1);
for i = 1:n+1
    h(i) = a1^n*a2^(i-1);
end

if nargout > 2
    B = zeros(n+1,1);
    for j = 1:n+1
        if mod(N,2)==0
            %B(j) = nchoosek(n,j-1)*((1/12)^(j-1)-2*(1/12+1/(24*N^2))^(j-1));
            B(j) = nchoosek(n,j-1)*((1/12)^(j-1)-2*(1/12+1/(24*N^2))^(j-1)...
                +0.25^(j-1)/N-(0.25-1/(6*s)-1/(6*N))^(j-1)/N);
        else
            %B(j) = nchoosek(n,j-1)*((1/12)^(j-1)-2*(1/12-1/(12*N^2))^(j-1));
            B(j) = nchoosek(n,j-1)*((1/12)^(j-1)-2*(1/12+1/(24*N^2))^(j-1)...
                +(0.25-1/(4*N^2))^(j-1)/N-(0.25-1/(6*s)-1/(6*N)-1/(4*N^2))^(j-1)/N);
        end
    end
    a1hat = a1-1;
    if mod(N,2)==0
        a2hat = (2*N-2)/(N*s-2);
    else
        a2hat = (2*N-2)/(N*s-2-3*s/N);
    end
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
[h,b,H,B] = CD_GWPoa_lhd(N,n,s);
z1 = norm(sum(H)-h);
z2 = norm(sum(B)-b);
%}