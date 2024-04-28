function [design,x] = Id2Design(q,x,mode)
% 此函数将各试验点对应字典排序的序号换为对应的设计，不计重复
% input:
%       q: s-by-1 向量，s个因子的水平数
%       x: n-by-1 binary vector, n 试验次数
%       mode: 1,输出按照字典排序，0，按照 x 的顺序排序(default)
% output:
%       design: 相应的设计，字典排序法, 水平从 0 开始

%narginchk(2, 2);%, nargin, 'struct'));
%nargoutchk(1, 1);%, nargout, 'struct'));

[s,m] = size(q);
if m~=1
    error('Id2Design:Input variable q must be a s-by-1 vector!\n');
end
if any( q~=round(q) ) || any( q <= zeros(size(q)) )
    error('Id2Design:The component of input variable must be positive integer!\n');
end
if any( x~=round(x) ) || any( x < zeros(size(x)) )
    error('Id2Design:The component of input variable must be positive integer or zero!\n');
end


n = length(x); % 试验次数
if nargin > 2 && mode ~= 0
    x = sort(x);
end

dd = ones(s,1);
for i=1:s-1
    dd(i)=prod(q(i+1:end));
end



design = zeros(n,s);
temp = x;
for k = 1:s
    design(:,k) = floor(temp/dd(k));
    temp = temp - design(:,k)*dd(k);
end

end