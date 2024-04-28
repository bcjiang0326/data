function [x,v] = Design2Bin(q,design)
% �˺��������ת��Ϊ��Ӧ���ֵ����򷨵� binary vector
% input:
%       q: s-by-1 ������s�����ӵ�ˮƽ��
%       design: ��Ӧ����ƣ�(�����ֵ����򷨣�, ˮƽ�� 0 ��ʼ
% output:
%       x: m-by-1 binary vector, m = prod(q) Ϊˮƽ�����
%       v: ����������Ӧ���ֵ����������
error(nargchk(2, 2, nargin, 'struct'));
error(nargoutchk(2, 2, nargout, 'struct'));

[s,m] = size(q);
if m~=1
    error('Design2Bin:Input variable q must be a s-by-1 vector!\n');
end
if any( q~=round(q) ) || any( q <= zeros(size(q)) )
    error('Design2Bin:The component of input variable must be positive integer!\n');
end
if size(design,2)~= s
    error('Design2Bin:The component of input variables design and q donot compatible!\n');
end 

m = prod(q);
x = zeros(m,1);
v = design(:,s);
temp = q(s);
for k = s-1:-1:1
    v = v + temp*(design(:,k));
    temp = temp*q(k);
end

x(v) = 1;
end