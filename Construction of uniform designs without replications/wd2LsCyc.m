function [x,bound,wd2] = wd2LsCyc(q,n,x)
% 搜索适用于Large-Scale问题的最优设计（WD2）,
% 外层迭代中令 x = [x(m);x(1:m-1)];
% Output:
%   x: solution
%   bound: 0.5*x'*H*x
%   wd2: the square of the wrap-around L2-discreption of x
% Input:
%   q: a vector which components are the numbers of levels for each factor
%   n: the runs of experiments
%   N: the number of repetitions
%   x: the initial solution.

% 2012/4/12 测试



error(nargchk(2, 4, nargin, 'struct'));
error(nargoutchk(0, 4, nargout, 'struct'));

[s,m] = size(q);
if m~=1
    error('wd2LsCyc:Input variable q must be a s-by-1 vector\n');
end

if any( q~=round(q) ) || any( q <= zeros(size(q)) )
    error('wd2LsCyc:The component of input variable must be positive integer\n');
end

if n~=round(n) || n <= 0
    error('wd2LsCyc:Input variable n must be a positive integer\n');
end

%{
if nargin >= 3
    if N~=round(N) || N <= 0
        error('wd2LsCyc:Input variable N must be a positive integer\n');
    end
else
    N = 1;
end
%}
    

global H; %H = getMatWD(q);
m = prod(q);
beq = mod(n,m);
k = (n-beq)/m;

if beq == 0
    x = k*ones(m,1);
    bound = 0.5*k*sum(sum(H));
    wd2 = -(4/3)^s + 2*bound/n^2;
    return;
end


if nargin < 3 || isempty(x)
    x = randperm(m)';
    v1 = x(1:beq);
    v0 = x(beq+1:end); 
    x(v1) = 1;
    x(v0) = 0;
else
    v1 = find(x==1); % father(v1(i)) == 1 
    v0 = find(x==0); % father(v0(i)) == 0
    if length(v1)+length(v0) ~= beq
        error('wd2LsCyc:The initial solutiaon must be m-by-1 binary vector\n');
    end
end

if length(v1)~=beq
    error('wd2LsCyc:The number of 1 in x0 must be beq\n');
end


epsilon = 1e-13;
bound = 0.5*sum(sum(H(v1,v1)));
v11 = v1;v00 = v0;

innerloop = 0;
while innerloop < max(q)
    innerloop = innerloop+1;
    % 迭代前产生一个初始的 x
    if innerloop == 2
        x(v11) = 1;
        x(v00) = 0;
    end
    if innerloop ~= 1
        x = [x(m);x(1:m-1)];
        v1 = find(x==1);
        v0 = find(x==0);
    end
    
    i = 0; j = 0; % father(i) = 1, father(j) = 0 , son(i) = 0, son(j) = 1 
    delta = -1;
    %iter = 0;
    while delta < -epsilon %&& iter < MaxIter 
        %iter = iter+1;
        delta = 0; % F(son) - F(father)
        for k1 = 1:length(v1)
            for l = 1:length(v0)
                if v1(k1)~=j && v0(l)~= i
                    if length(v1) ~= 1
                        temp_delta = -( 2*sum(H(v1(k1),v1)) - H(v1(k1),v1(k1)) ) + ...
                            ( 2*sum(H(v0(l),[v1(1:k1-1);v1(k1+1:end);v0(l)])) - H(v0(l),v0(l)) );
                    else
                        temp_delta = H(v0(l),v0(l)) - H(v1,v1);                        
                    end
                    if temp_delta < delta
                        delta = temp_delta;
                        k_ = k1;
                        l_ = l;
                    end
                end
            end
        end
        %if son is better then update i, j, v1, v0
        if delta < 0
            i = v1(k_); j = v0(l_);
            temp = v1(k_);
            v1(k_) = v0(l_);
            v0(l_) = temp;
        end
    end
    bound_ = 0.5*sum(sum(H(v1,v1))); %用 bound—— 计算此次算出的x对应的bound
    if bound_ < bound - epsilon
        bound = bound_; % 用 bound 记录到目前为止最小的 bound
        v11 = v1;   % 用v11记录到目前为止最优的 x 中取 1 的位置
        v00 = v0;   % 用v00记录到目前为止最优的 x 中取 1 的位置
    end
end
        
x(v11) = 1;
x(v00) = 0;
if n > m         
    x = x + ones(m,1)*k;
end

wd2 = -(4/3)^s + x'*H*x/n^2;
end
