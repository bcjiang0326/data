function design = Bin2Design(q,x)
% �˺������ֵ����򷨻�õ� x ֵת��Ϊ��Ӧ�����
% input:
%       q: s-by-1 ������s�����ӵ�ˮƽ��
%       x: m-by-1 binary vector, m = prod(q) Ϊˮƽ�����
% output:
%       design: ��Ӧ����ƣ��ֵ����򷨣���ƴ� 0 ��ʼ

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

n = sum(x); % �������
temp = find(x); % ȷ���������������ֵ������е����

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