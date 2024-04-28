function y = LB_aveDisc_h(c,M,m)
%20190416 by Bochuan Jiang
%Ŀ�꺯��: sum_{i=1}^m c^{x_i}�����У�c > 1���� 
%Լ��: (1)sum_{i=1}^m x_i = M
%      (2)x_i �Ǹ�����
%�����򷵻�����Ŀ�꺯��
%INPUT:
%   c: c>1
%   M: ������
%   m: ������
%OUTPUT:
%   y: ����Ŀ�꺯��
if c <= 1
    error('c > 1!\n');
end
if M~=floor(M) || M<=0 || m~=floor(m) || m<=0
    error('Wrong M or m!\n');
end
a = floor(M/m);
k1 = m*a+m-M;
k2 = m-k1;
y = c^a*k1+c^(a+1)*k2;
end