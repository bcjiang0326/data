function y = WD2_From_comp(q,n,wd)
% ͨ�� D �Ĳ���� bar_D �� WD^2������ WD^2(D)
% Input: 
%       q: m-by-1 vector ��Ԫ��Ϊ�����ӵ�ˮƽ��
%       n: D ������������� N-n Ϊ bar_D ���������
%       wd: WD^2(bar_D)
% Output:
%       y: = WD^2(D)
m = length(q);
N = prod(q);

a1 = -(4/3)^m;
a2 = (2*n-N)/n^2;
for k = 1:m
    a2 = a2*(4*q(k)/3+1/(6*q(k)));
end
a3 = (N-n)^2/n^2*(wd+(4/3)^m);
y = a1+a2+a3;

end


