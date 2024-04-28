function [c0,c1,c2] = aveDisc_c1c2(s,flag)
%20190416 by Bochuan Jiang
%在仅考虑水平置换情形下，计算average discrepancy 中的 c1(s) 和 c2(s)
%Input:
%   s: s > 1, 正整数
%   flag: 'CD'(default),'WD','MD'
%Output:
%   c0: 一个常数 c0 = int_0^1 int_0^1 f(x,y) dx dy
%   c1: c1(s)
%   c2: c2(s)
if any(s<=1) || any(s~=floor(s))
    error('aveDisc_c1c1: Wrong s!\n');
end
if strcmp(flag,'CD')
    c0 = 13/12;
elseif strcmp(flag,'WD')
    c0 = 4.0/3.0;
elseif strcmp(flag,'MD')
    c0 = 19/12;
else
    error('aveDisc_c1c2: Wrong flag!\n');
end
c1 = zeros(size(s));
for i = 1:length(s)
    for z1 = 0:s(i)-2
        for z2 = z1+1:s(i)-1
            a = (z1+0.5)/s(i);
            b = (z2+0.5)/s(i);
            c1(i) = c1(i)+kernel_f(a,b,flag);
        end
    end
    c1(i) = c1(i)*2;
end
c2 = zeros(size(s));
for i = 1:length(s)
    for z = 0:s(i)-1
        a = (z+0.5)/s(i);
        c2(i) = c2(i)+kernel_f(a,a,flag);
    end
    c2(i) = c2(i)*(s(i)-1)/c1(i);
end
end

