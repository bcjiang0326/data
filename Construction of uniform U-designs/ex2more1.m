% 此脚本目的：计算 OA(36,2^3 3^3,2) 的 LHD
epsilon = 1e-10; T0 = 1e-2;  T1 = 1e-6;
InIter = 1e4; OutIter = 100*InIter; Reps = 5;

OA = importdata('OA from web/MA.36.3.12.2.11.txt');
[N,ncols] = size(OA);
q = zeros(ncols,1);
for i = 1:ncols
    q(i) = length(unique(OA(:,i)));
end
[q,id] = sort(q);
OA = OA(:,id);

q_oa = [2;2;2;3;3;3];
id = [1 2 4 17 18 19];
D = OA(:,id);
n10 = sum(q_oa==2);
n20 = sum(q_oa==3);
n = length(q_oa);
[~, Disc_vec0, L0] = TA_MD_LHD(D,q_oa,OutIter,InIter,T0,T1,Reps,0,1);
%Lname = strcat('MA LHD 20151024 MD oa36 mwa ta/L_mwa',int2str(n10),int2str(n20),'_',int2str(4));
out= fopen(Lname,'w');
for i = 1:N
    fprintf(out,'%d ',L0(i,:));
    fprintf(out,'\n');
end
fclose(out);

ex2more
