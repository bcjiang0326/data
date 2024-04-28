oa = importdata('../Web OA/MA.72.4.1.6.5.3.8.2.7.txt');
allDisc = {'WD';'CD';'MD'};

N = size(oa,1);
Part1 = [oa(:,1),rand_U_Type_design(4*ones(6,1),N)];
Part2 = [oa(:,2:6),rand_U_Type_design(6*ones(2,1),N)];
Part3 = oa(:,7:14);
oa0 = [Part1,Part2,Part3];
n = size(oa0,2);
q0 = zeros(n,1);
for j = 1:n
    q0(j) = length(unique(oa0(:,j)));
end
[isgrouped,uq] = isgrouped_ByLevNum(q0);
n1 = sum(q0==uq(1));
n2 = sum(q0==uq(2));
n3 = sum(q0==uq(3));



Jset = [2:6,13,14]';
InIter1 = 1000;
AllIter1 = 1000*InIter1;
outfile1 = 'MA.72.Supersaturated.txt';
InIter2 = 1000;
AllIter2 = 100*InIter1;
outfile2 = 'T_Disc_Disc.txt';

T1 = 0; t = cputime;
[aveWD1, aveWD1vec, D0, isLB] = TA_MAD_Local(oa0,Jset,AllIter1,InIter1,1e-3,1e-6,'WD',1,outfile1);
[WD1, WD1vec, D1] = TA_uniform_FF_design_WD(D0,AllIter2,InIter2,1e-3,1e-6,1,1);
T1 = T1+cputime-t;
fprintf('%,8f\n',WD1);

figure
subplot(1,2,1);
data1 = importdata(outfile1);
plot(data1);
xlabel('Iterations');
ylabel('Average discrepancy');
title('Stage 1');
subplot(1,2,2);
data2 = importdata(outfile2);
plot(data2(:,3));
xlabel('Iterations');
ylabel('Discrepancy');
title('Stage 2');







