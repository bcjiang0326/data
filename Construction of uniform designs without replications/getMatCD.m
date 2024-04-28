function [A f] = getMatCD(q)
% CD^2(U) = (13/12)^s-2/n*f'*y + 1/n^2*y'*A*y
% Input:
% q: s-by-1 vector with positive integer component
%
% Output:
% A: A = A_1#...#A_s where # means Kronecker tensor product and A_k is
%    q_k-by-q_k matrix with A_k(i,j)=1+0.25*|(2*i+1-q_k)/q_k|
%    +0.25*|(2*j+1+q_k)/q_k|-0.5*|(i-j)/q_k|
% f: f = f1#...#f_s where fk(i) = 1+0.25*|(2*i+1-q_k)/q_k|
%        -0.125|(2*i+1-q_k)/q_k|^2

[n,m] = size(q);
if min(m,n)~=1
    error('getMatCD:Input variable must be a vector\n');
end

if any(abs(q-round(q))) || any(abs(q-abs(q))) || any( q == zeros(size(q)) )
    error('getMatCD:The component of input variable must be positive integer\n');
end

n = max(n,m); A = 1;
f = 1;
for s = 1:n
    As = zeros(q(s),q(s));
    fs = zeros(q(s),1);
    for i = 1:q(s)
        for j = 1:i
            As(i,j)=1+0.25*abs( (2*i-1-q(s))/q(s) )+...
                0.25*abs( (2*j-1-q(s))/q(s) )-...
                0.5*abs((i-j)/q(s));
            if j~=i 
                As(j,i) = As(i,j);
            end
        end
        fs(i) =  1+0.25*abs((2*i-1-q(s))/q(s))...
            -0.125*abs((2*i-1-q(s))/q(s))^2;
    end
    A = kron(A,As);
    f = kron(f,fs);
end
end
            