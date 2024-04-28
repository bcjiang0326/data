function [D0,Disc0,id0] = Disc_sub(D,q0,flag)
%20181208 by Bochuan Jiang
%按照 average discrepancy under level-permutating 最小化原则, 采用穷举法
%从 D 中选取子设计D0, 而子设计的各列水平数为q0(i).
%INPUT:
%   D: 一个非对称 OA
%   q0: n-by-1 向量，记录OA各列的水平数
%   flag: 'CD'(default),'WD','MD'
%OUTPUT:
%   D0: 得到的设计
%   Disc: 得到的设计的 discrepancy
%   id0: D(:,id0) = D0;

epsilon = 1e-10;
if nargin < 3
    flag = 'CD';
end
if strcmp(flag,'CD')
    Fun = @CD2_value;
elseif strcmp(flag,'WD')
    Fun = @WD2_value;
elseif strcmp(flag,'MD')
    Fun = @MD2_value;
else
    error('Wrong flag\n');
end

n = size(D,2);
q = zeros(n,1);
for i = 1:n
    q(i) = length(unique(D(:,i)));
end
[isgrouped,uq,m] = isgrouped_ByLevNum(q);
if ~isgrouped
    error('q is not grouped by number of levels!\n');
end

m0 = zeros(length(m),1);
q0_ = zeros(0,1);
for k = 1:length(uq)
    m0(k) = sum(q0==uq(k));
    q0_  = cat(1,q0_,uq(k)*ones(m0(k),1));
end
if sum(m0)~=length(q0)
    error('q0 and D are not match!\n');
end
q0 = q0_;
u = m0~=0; % m0 不为 0 的下角标

dd = zeros(1,m0(1));
for i = 2:length(m0)
    dd = cat(2,dd,(dd(end)+m(i-1))*ones(1,m0(i)));
end
id = 1:m0(1);
for i = 2:length(m)
    id = cat(2,id,(1:m0(i)));
end

Disc0 = inf;
while id(end)~=-1     
    Disc = Fun(D(:,id+dd),q0);
    %{
    fprintf('%.4e ',aveCD); fprintf('%d ',id+dd); fprintf('\n');
    fprintf('%.4e ',A(3:end)); fprintf('\n');
    fprintf('\n');
    %}
    if Disc < Disc0-epsilon
        Disc0 = Disc;
        id0 = id;
    end
    id = choosenextgroup(m(u),m0(u),id);
end
id0 = id0+dd;
D0 = D(:,id0);

end

%{
oa = importdata('OA from web/MA.36.3.7.6.3.finney.txt');
q0 = [3;3;3;3;6];
[D0,Disc0,id0] = Disc_sub(oa,q0);
%}

