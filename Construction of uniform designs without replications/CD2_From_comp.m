function y = CD2_From_comp(q,id_bar_D)
% ͨ�� D �Ĳ���� bar_D �� CD^2������ CD^2(D)��
% n Ϊ D ������������� N-n Ϊ bar_D �����������
% Input: 
%       q: m-by-1 vector ��Ԫ��Ϊ�����ӵ�ˮƽ��
%       id_bar_D: bar_D �� rank vector
% Output:
%       y: = CD^2(D)
[~,D] = compID(id_bar_D,q);
y = CD2_value( D, q );
end


