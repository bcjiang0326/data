function [D,v] = sortDesign(D,q,mode)
% 20121207
% 将设计按照字典排序法排序
% INPUT:
%       D: 设计，水平从0开始
%       q: 各因子水平数
%       mode: 排序方式
%           'ascend': 升序。默认
%           'descend': 降序
% OUTPUT:
%       D: 排序后的设计
%       v: 字典排序的序号数 D_in(v)=D_out
v = Design2Id(q,D);
if nargin > 2
    [v,ix] = sort(v,mode);
else
    [v,ix] = sort(v);
end
D = D(ix,:);
