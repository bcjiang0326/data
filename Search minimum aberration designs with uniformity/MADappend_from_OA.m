function [D0,aveDisc0,id0] = MADappend_from_OA(D,q0,flag,filename)
%20190802 by Bochuan Jiang
%按照 average discrepancy under level-permutating 最小化原则, 采用穷举法
%从 D 中选取子设计D0, 使得列并置[D,D0]的average discrepancy 最小。
%子设计的各列水平数为记为q0(i).
%INPUT:
%   D: 一个非对称 OA
%   q0: n-by-1 向量，记录OA各列的水平数
%   flag: 'CD'(default),'WD' or 'MD'
%   filename: 提供文件名的话，则将结果写入文件
%OUTPUT:
%   D0: 得到的设计
%   aveDisc: 得到的设计的 average discrepancy
%   id0: D(:,id0) = D0;

epsilon = 1e-10;
if isempty(q0)
    error('Wrong q0!\n');
end
if nargin < 3
    flag = 'CD';
end
if ~strcmp(flag,'CD') && ~strcmp(flag,'WD') && ~strcmp(flag,'MD')
    error('Wrong flag!\n');
end
if nargin < 4
    iswrite = 0;
else
    iswrite = 1;
    outfile = fopen(filename,'w');
end

n = size(D,2);
q = zeros(n,1);
for i = 1:n
    q(i) = length(unique(D(:,i)));
end
[isgrouped,uq,m] = isgrouped_ByLevNum(q);
if ~isgrouped
    error('MADsub_from_OA: q is not grouped by number of levels!\n');
end

m0 = zeros(length(m),1);
q0_ = zeros(0,1);
for k = 1:length(uq)
    m0(k) = sum(q0==uq(k));
    q0_  = cat(1,q0_,uq(k)*ones(m0(k),1));
end
if sum(m0)~=length(q0)
    error('MADsub_from_OA: q0 and D are not match!\n');
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

aveDisc0 = inf;
while id(end)~=-1
    q_ = [q;q0];
    D_ = [D,D(:,id+dd)];
    if ~issorted(q_)
        [q_,a] = sort(q_);
        D_ = D_(:,a);
    end
    aveDisc = aveDisc_LevelPerm(D_,q_,flag);
    if iswrite
        fprintf(outfile,'%.8e ',aveDisc);
        fprintf(outfile,'%d ',id+dd);
        fprintf(outfile,'\n');
    end
    if aveDisc < aveDisc0-epsilon
        aveDisc0 = aveDisc;
        id0 = id;
    end
    id = choosenextgroup(m(u),m0(u),id);
end
id0 = id0+dd;
D0 = D(:,id0);

if iswrite
    fclose(outfile);
end

end

%{
oa = importdata('OA from web/MA.36.3.7.6.3.finney.txt');
q0 = [3;6];
[D0,aveDisc0,id0] = MADappend_from_OA(oa,q0);
%}

%{
oa = importdata('OA from web/oa.27.13.3.2.txt');
q0 = [3;3];
[D0,aveDisc0,id0] = MADappend_from_OA(oa,q0,'WD');
fprintf('%.6f\n',aveDisc0);
%}

