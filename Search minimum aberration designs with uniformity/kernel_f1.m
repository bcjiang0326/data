function z = kernel_f1(x,flag)
%20190416 by Bochuan Jiang
%计算 discrepancy 对应的 reproducing kernel 中的 f1(x)=inf_0^1 f(x,y) dy
%Input:
%   x: n-by-1 vector, 每个元素属于[0,1]
%   flag: 'CD'(default),'WD','MD'
%Output:
%   z: f1(x)
if any(x<0) || any(x>1) || size(x,2)~=1 
    error('kernel_f1: wrong x !\n');
end

if nargin < 2
    flag = 'CD';
end
if strcmp(flag,'CD')
    z = 1+(abs(x-0.5)-(x-0.5).^2)/2; 
elseif strcmp(flag,'WD')
    z = 4.0/3.0*ones(size(x));
elseif strcmp(flag,'MD')
    z = 5.0/3.0-(abs(x-0.5)+(x-0.5).^2)/4;
else
    error('kernel_f1:Wrong flag!\n');
end
end

