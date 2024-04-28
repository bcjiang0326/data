function [design,x] = Id2Design(q,x,mode)
% �˺�������������Ӧ�ֵ��������Ż�Ϊ��Ӧ����ƣ������ظ�
% input:
%       q: s-by-1 ������s�����ӵ�ˮƽ��
%       x: n-by-1 binary vector, n �������
%       mode: 1,��������ֵ�����0������ x ��˳������(default)
% output:
%       design: ��Ӧ����ƣ��ֵ�����, ˮƽ�� 0 ��ʼ

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


n = length(x); % �������
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