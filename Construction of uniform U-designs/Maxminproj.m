function y = Maxminproj(L,dims,p)
% 20151019 Xiaopang
% �����ض�ά��ͶӰ���ƽ�� Maxminֵ
% L ������һ��Ԫ��Ϊ{0,...,N-1} ������������
if nargin < 3
    p = 3;
end
if size(dims,1)~=1
    error('dims must by a row vector!\n');
end

[N,n] = size(L);
L = (L+0.5)/N;

y = zeros(size(dims));
for i = 1:length(dims)
    d = dims(i);
    y1 = 0;
    id = 1:d;
    while id(1)~=-1
        y1 = y1+MaxminValue(L(:,id),p);
        id = nchoosek_next(id,n,d);
    end
    y(i) = y1/nchoosek(n,d);
end
end