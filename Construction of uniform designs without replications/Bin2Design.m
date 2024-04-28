function design = Bin2Design(q,x)
% 此函数将字典排序法获得的 x 值转换为对应的设计
% input:
%       q: s-by-1 向量，s个因子的水平数
%       x: m-by-1 binary vector, m = prod(q) 为水平组合数
% output:
%       design: 相应的设计，字典排序法，设计从 0 开始

%narginchk(2, 2);%, nargin, 'struct'));
%nargoutchk(1, 1);%, nargout, 'struct'));

[s,m] = size(q);
if m~=1
    error('Bin2Design:Input variable q must be a s-by-1 vector!\n');
end
if any( q~=round(q) ) || any( q <= zeros(size(q)) )
    error('Bin2Design:The component of input variable must be positive integer!\n');
end
if length(x(x==1)) + length(x(x==0)) ~= length(x) 
    error('Bin2Design:The input variable x must be m-by-1 binaries!\n');
end

n = sum(x); % 试验次数
temp = find(x); % 确定各个试验在总字典排序中的序号

dd = ones(s,1);
for i=1:s-1
    dd(i)=prod(q(i+1:end));
end
design = zeros(n,s);
for k = 1:s
    design(:,k) = floor(temp/dd(k));
    temp = temp - design(:,k)*dd(k);
end

end