function [D0,id0,aveDisc0] = MADsub_from_OA(D,q0,flag,filename)
%20181208 by Bochuan Jiang
%���� average discrepancy under level-permutating ��С��ԭ��, ������ٷ�
%�� D ��ѡȡ�����D0, ������Ƶĸ���ˮƽ��Ϊq0(i).
%INPUT:
%   D: һ���ǶԳ� OA
%   q0: n-by-1 ��������¼OA���е�ˮƽ��
%   flag: 'CD'(default),'WD' or 'MD'
%   filename: �ṩ�ļ����Ļ����򽫽��д���ļ�
%OUTPUT:
%   D0: �õ������
%   id0: D(:,id0) = D0
%   aveDisc: �õ�����Ƶ� average discrepancy


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
u = m0~=0; % m0 ��Ϊ 0 ���½Ǳ�

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
    aveDisc = aveDisc_LevelPerm(D(:,id+dd),q0,flag);
    if iswrite
        fprintf(outfile,'%.8e ',aveDisc);
        fprintf(outfile,'%d ',id+dd);
        fprintf(outfile,'\n');
    end
    %{
    fprintf('%.4e ',aveDisc); fprintf('%d ',id+dd); fprintf('\n');
    %}
    if aveDisc < aveDisc0-epsilon
        aveDisc0 = aveDisc;
        id0 = id;
        %fprintf('%d ',id0); fprintf('\n');
        %fprintf('%.8f,',aveDisc0); fprintf('%d ',id0+dd); fprintf('\n');
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
oa = importdata('../Web OA/MA.36.3.7.6.3.finney.txt');
q0 = [3;3;3;3;6];
[D0,id0,aveDisc0] = MADsub_from_OA(oa,q0,'CD','out.txt');
%}


%{
oa = importdata('../Web OA/MA.36.3.12.2.11.txt');
q0 = [3;3;3;2;2;2];
[D0,id0,aveDisc0] = MADsub_from_OA(oa,q0,'CD','out.txt');
%}

