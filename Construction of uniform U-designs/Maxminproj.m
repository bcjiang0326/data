function y = Maxminproj(L,dims,p)
% 20151019 Xiaopang
% 搜索特定维度投影设计平均 Maxmin值
% L 必须是一个元素为{0,...,N-1} 的拉丁超立方
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