function [D,y,Disc] = DesignOnWeb(levels,s,n,flag)
% 20130127
% 获取网站设计
% INPUT:
%       levels: 水平数
%       s: 因子数
%       n: 试验数
%       flag: 0(default),CD; 1, WD
% OUTPUT:
%       D: 设计，按字典排序
%       y: rank vector
%       Disc: square of discrepancy

if nargin < 4
    flag = 0;
end

if flag
    fname='WD2';
else
    fname='CD2';
end

if levels == 3
        fname = strcat(fname,'/Level 3/',int2str(s),'_',int2str(n),'.txt');
    elseif levels == 4
        fname = strcat(fname,'/Level 4/',int2str(s),'_',int2str(n),'.txt');
    elseif levels == 5
        fname = strcat(fname,'/Level 5/',int2str(s),'_',int2str(n),'.txt');
    elseif levels == 6
        fname = strcat(fname,'/Level 6/',int2str(s),'_',int2str(n),'.txt');
end
D = importdata(fname)-1;
q = levels*ones(s,1);
[D,y] = sortDesign(D,q);

if nargout > 2
    if flag
        Disc = WD2_value(D,q);
    else
        Disc = CD2_value(D,q);
    end
end

end