function [y,y1,y2,y2d] = MD2_value( U, q )
% 2013 0403
% q 水平因子的各水平为: 0,...,q-1
% INPUT:
%       U: 设计
%       q: 各因子的水平
% OUTPUT:
%       y: Squared discrepancy
%       y1: 一次项求和的绝对值
%       y2: 二次项 i不等于j 部分 求和
%       y2d: 二次项 i等于j 部分 求和 

[n,s] = size(U);

if length(q) ~= s || size(q,2)~= 1
    error('q and U do not match each other!\n');
end
q = q';

y2 = 0;
for i = 1:n-1
    for j = i+1:n
        y2 = y2 + prod( 1.875 - 0.125*abs(2*U(i,:)+1-q)./q...
            - 0.125*abs(2*U(j,:)+1-q)./q...
            - 0.75*abs(U(i,:)-U(j,:))./q...
            + 0.5*((U(i,:)-U(j,:))./q).^2);
    end
end
y2 = y2*2/n^2;

y2d = 0;
for i = 1:n
    y2d = y2d + prod( 1.875 - 0.25*abs(2*U(i,:)+1-q)./q );
end
y2d = y2d/n^2;

y1 = 0;
for i = 1:n
    y1 = y1 + prod( 5/3 - 0.125*abs(2*U(i,:)+1-q)./q...
            - 0.0625*((2*U(i,:)+1-q)./q).^2 );
end
y1 = y1*2/n;

y = (19/12)^s - y1 + y2 + y2d;

end
%{
% 测试代码
D = [0,1,2,3,4,5,6; 4,1,6,3,0,5,2]';
q = [7;7];
[y,y1,y2,y2d] = MD2_value( D, q );
%}