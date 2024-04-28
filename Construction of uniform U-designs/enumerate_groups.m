%20160512 by xiaopang
%���ض������ɸ������г�����ÿ�����ϳ�ȡȷ�����������������ȫ���ĳ�������
%INPUT:
%   M: ����Ԫ��Ϊ�������ϵ�Ԫ�ظ���
%   M0: ����Ԫ��Ϊ�������ϵĴ���ȡ��Ԫ�ظ���
%OUTPUT:
%   vec: ÿ��Ϊһ����������
function vec = enumerate_groups(M,M0)
N = length(M);
id = zeros(1,0);
for i = 1:N
    id = cat(2,id,1:M0(i));
end
vec = zeros(0,length(id));
while id(end)~=-1
    vec = cat(1,vec,id);
    id = choosenextgroup(M,M0,id);
end
end

%{
%���Գ���
M = [3;4;5];
M0 = [1;3;4];
vec = enumerate_groups(M,M0);

%}