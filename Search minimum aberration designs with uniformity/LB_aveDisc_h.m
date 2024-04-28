function y = LB_aveDisc_h(c,M,m)
%20190416 by Bochuan Jiang
%目标函数: sum_{i=1}^m c^{x_i}，其中（c > 1）。 
%约束: (1)sum_{i=1}^m x_i = M
%      (2)x_i 非负整数
%本程序返回最优目标函数
%INPUT:
%   c: c>1
%   M: 正整数
%   m: 正整数
%OUTPUT:
%   y: 最优目标函数
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