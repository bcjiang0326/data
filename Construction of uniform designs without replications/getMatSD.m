function [A f] = getMatSD(q)
% SD^2(U) = (4/3)^s-2/n*f'y + 1/n^2*y'*A*y
% Input:
% q: n-by-1 vector with positive integer component
%
% Output:
% A: A = A_1#...#A_s where # means Kronecker tensor product and A_k is
%    q_k-by-q_k matrix with A_k(i,j)=2-2|i-j|/q_k
% f: f = f1#...#f_s where f_k(i) = 1+(2i-1)/q_k-(2i-1)^2/(2q_k^2)

[n,m] = size(q);
if min(m,n)~=1
    error('GetA:Input variable must be a vector\n');
end

if any(abs(q-round(q))) || any(abs(q-abs(q))) || any( q == zeros(size(q)) )
    error('GetA:The component of input variable must be positive integer\n');
end

n = max(n,m); A = 1;
f = 1;
for s = 1:n
    As = zeros(q(s),q(s));
    fs = zeros(q(s),1);
    for i = 1:q(s)
        for j = 1:i
            As(i,j) = 2-2*abs(i-j)/q(s);
            if(j~=i)
                As(j,i) = As(i,j);
            end
        end
        fs(i) =  1 + (2*i-1)/q(s)-(2*i-1)^2/(2*q(s)^2);
    end
    A = kron(A,As);
    f = kron(f,fs);
end
end
            