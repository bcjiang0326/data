function [ID1,D1] = compID(ID0,q)
% 20121211
% �о�����Ƶ��ֵ��������ֵ
% ע��ԭ��Ʊ������ظ�
% INPUT:
%       ID0: ԭ��Ƶ��ֵ��������ֵ���������ظ�
%       q: s-by-1 vector, ��������ˮƽ��
% OUTPUT:
%       ID1: ����Ƶ��ֵ��������ֵ

if size(q,2) ~= 1 || size(ID0,2)~= 1
    error('compID: q and ID0 must be vectors!\n');
end

ID0 = sort(ID0);
N = prod(q);
n = length(ID0);
ID1 = zeros(N-n,1);
k = 1; begin = 0;
for i = 1:n
    for j = begin:ID0(i)-1
        ID1(k) = j;
        k = k+1;
    end
    begin = ID0(i)+1;
end

for j = begin:(N-1)
    ID1(k) = j;
    k = k+1;
end

if nargout > 1
    D1 = Id2Design(q,ID1);
end
end

