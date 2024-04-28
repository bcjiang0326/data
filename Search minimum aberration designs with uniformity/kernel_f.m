function z = kernel_f(x,y,flag)
%20190416 by Bochuan Jiang
%计算 discrepancy 对应的 reproducing kernel 中的 f(x,y)
%Input:
%   x: n-by-1 vector, 每个元素属于[0,1]
%   y: n-by-1 vector, 每个元素属于[0,1]
%   flag: 'CD'(default),'WD','MD'
%Output:
%   z: f(x,y)
if any(x<0) || any(x>1) || any(y<0) || any(y>1) || size(x,2)~=1 || size(y,2)~=1
    error('kernel_f: wrong x or y!\n');
end
if ~strcmp(flag,'CD') && ~strcmp(flag,'WD') && ~strcmp(flag,'MD')
    error('kernel_f: wrong flag!\n');
end
if nargin < 2
    flag = 'CD';
end
if strcmp(flag,'CD')
    z = 1+(abs(x-0.5)+abs(y-0.5)-abs(x-y))/2; 
elseif strcmp(flag,'WD')
    z = 1.5-abs(x-y)+(x-y).^2;
elseif strcmp(flag,'MD')
    z = 1.875-(abs(x-0.5)+abs(y-0.5)+3*abs(x-y)-2*(x-y).^2)/4;
else
    error('kernel_f:!!!!\n');
end
end

