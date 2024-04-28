function [D0,id0,A0] = GMAsub_from_OA(D,q0,filename)
%20181208 by Bochuan Jiang
%按照 GMA 准则, 采用穷举法从混水平设计 D 中选取子设计D0, 而子设计的各列水平数为q0(i).
%INPUT:
%   D: 一个非对称 OA
%   q0: n-by-1 向量，记录OA各列的水平数
%   flag: 'CD'(default),'WD' or 'MD'
%   filename: 提供文件名的话，则将结果写入文件
%OUTPUT:
%   D0: 得到的设计
%   id0: D(:,id0) = D0;
%   A0: 得到的设计对应的GWP

epsilon = 1e-10;
if isempty(q0)
    error('Wrong q0!\n');
end
if nargin < 3
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

A0 = inf(1,length(q0));
while id(end)~=-1
    if all(q0==q0(1))
        A = GWP(D(:,id+dd));
    else
        A = GWP_asym(D(:,id+dd));
    end
    if iswrite
        %fprintf(outfile,'%.8e ',aveDisc);
        fprintf(outfile,'%d ',id+dd);
        fprintf(outfile,'%.12e ',A);
        fprintf(outfile,'\n');
    end
    %{
    fprintf('%d ',id+dd); fprintf('\n');
    fprintf('%.4e ',A(3:end)); fprintf('\n');
    fprintf('\n');
    %}
    for j = 1:length(q0)
        if abs(A(j) - A0(j))>epsilon
            if A(j) < A0(j)-epsilon
                A0 = A;
                id0 = id;
                %fprintf('%.4f ',A0(3:end));fprintf('\n');
            end
            break;
        end
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
[D0,id0,A0] = GMAsub_from_OA(oa,q0,'out.txt');
%}

%{
oa = importdata('../Web OA/MA.36.2.11.3.12.txt');
q0 = [2;2;2;3;3;3;];
[D0,id0,A0] = GMAsub_from_OA(oa,q0,'out.txt');
%}

