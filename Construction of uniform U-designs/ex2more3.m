% 本程序去掉 MWA 算法中的 n10 = 3 n20 = 4 的最小MD对应的设计，换入具有更小MD值的设计
% 为保持平均 MD 保持不变，还要调整另外的设计
ex2more

OA = importdata('OA from web/MA.36.3.12.2.11.txt');
[N,ncols] = size(OA);
q = zeros(ncols,1);
for i = 1:ncols
    q(i) = length(unique(OA(:,i)));
end
[q,id] = sort(q);
OA = OA(:,id);

epsilon = 1e-10; T0 = 1e-2;  T1 = 1e-6;
InIter = 1e4; OutIter = 100*InIter; Reps = 1000;

%{
q_oa =[2;2;2;3;3;3;3];
id=[1 2 3 16 21 22 23];
D = OA(:,id);
n10 = sum(q_oa==2);
n20 = sum(q_oa==3);
n = length(q_oa);
q_lh = N*ones(n,1);
Disc2 = zeros(Reps,1);
L2 = cell(Reps,1);
for rep = 1:Reps
    [Disc2(rep), ~, L2{rep}] = TA_MD_LHD(D,q_oa,OutIter,InIter,T0,T1,1,0,1);
end
%}

Disc2 = importdata('Disc2.mat');
L2 = importdata('L2.mat');

[~,id0] = sort(Disc0);
aveDisc0 = mean(Disc0);
[Disc2,a] = sort(Disc2);
L2 = L2(a);
a = Disc2>0;
Disc2 = Disc2(a);
L2 = L2(a);

i0 = 0;
if Disc2(1)<min(Disc0)-epsilon
    i0 = 1;
end


if i0 > 0
    diff = inf;
    j0 = 0;
    for j = i0+1:Reps
        diff_ = abs(mean([Disc2(i0);Disc0(id0(2:end-1));Disc2(j)])-aveDisc0);
        if diff_ < diff-epsilon
            diff = diff_;
            j0 = j;
        end
    end
end




