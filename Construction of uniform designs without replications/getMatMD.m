function [A f] = getMatMD(q)
% MD^2(U) = (19/12)^s-2/n*f'*y + 1/n^2*y'*A*y
% Input:
% q: s-by-1 vector with positive integer component
%
% Output:
% A: A = A_1#...#A_s where # means Kronecker tensor product and A_k is
%        a q_k-by-q_k matrix 
% f: f = f1#...#f_s 

[n,m] = size(q);
if min(m,n)~=1
    error('getMatMD:Input variable must be a vector\n');
end

if any(abs(q-round(q))) || any(abs(q-abs(q))) || any( q == zeros(size(q)) )
    error('getMatMD:The component of input variable must be positive integer\n');
end

n = max(n,m); A = 1;
f = 1;
for s = 1:n
    As = zeros(q(s),q(s));
    fs = zeros(q(s),1);
    for i = 1:q(s)
        for j = 1:i
            As(i,j) = 1.875 - ( 0.125*abs( (2*i-1-q(s)) ) ...
                + 0.125*abs( (2*j-1-q(s)) ) ...
                + 0.75*abs( (i-j) ) ...
                - 0.5/q(s)*(i-j)^2 )/q(s);
            if j~=i 
                As(j,i) = As(i,j);
            end
        end
        fs(i) =  5/3 - ( 0.125*abs( (2*i-1-q(s)) ) ...
            + 0.0625*(2*i-1-q(s))^2/q(s) ) / q(s);
    end
    A = kron(A,As);
    f = kron(f,fs);
end
end
%{
% 测试程序
D = [0,1,2,3,4,5,6; 4,1,6,3,0,5,2]';
q = [7;7]; 
[A f] = getMatMD(q);
N = prod(q); [n,s] = size(D); y = zeros(N,1);
ID = Design2Id(q,D);
y(ID+1) = 1;
MD2 = (19/12)^s-2/n*f'*y + 1/n^2*y'*A*y

%{
% 此代码解释 MD 下补设计规律
q = randi(10);
[A,b] = getMatMD(q);
y = A*ones(size(b))-q*b-7/48/q;
if mod(q,2)~=0
    y = y-1/16/q;
end
%}
%}

