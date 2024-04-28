function [D,v] = sortDesign(D,q,mode)
% 20121207
% ����ư����ֵ���������
% INPUT:
%       D: ��ƣ�ˮƽ��0��ʼ
%       q: ������ˮƽ��
%       mode: ����ʽ
%           'ascend': ����Ĭ��
%           'descend': ����
% OUTPUT:
%       D: ���������
%       v: �ֵ����������� D_in(v)=D_out
v = Design2Id(q,D);
if nargin > 2
    [v,ix] = sort(v,mode);
else
    [v,ix] = sort(v);
end
D = D(ix,:);
