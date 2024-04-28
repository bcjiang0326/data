function [A f] = getMatWD(q)
% Input:
% q: n-by-1 vector with positive integer component
%
% Output:
% A: A = A_1#...#A_s where # means Kronecker tensor product and A_i means
%    q_i-by-q_i matrix with A_i(i,j) = 1.5 - |i-j|(q_i-|i-j|)/q_i^2

[n,m] = size(q);
if min(m,n)~=1
    error('GetA:Input variable must be a vector\n');
end

if any(abs(q-round(q))) || any(abs(q-abs(q))) || any( q == zeros(size(q)) )
    error('GetA:The component of input variable must be positive integer\n');
end

n = max(n,m); A = 1;
for s = 1:n
    As = 1.5*diag(ones(q(s),1));
    for i = 2:q(s)
        for j = 1:i-1
            As(i,j) = 1.5 - (i-j)*(q(s)-(i-j))/q(s)^2;
            As(j,i) = As(i,j);
        end
    end
    A = kron(A,As);
end
f = zeros(length(A),1);
end
            