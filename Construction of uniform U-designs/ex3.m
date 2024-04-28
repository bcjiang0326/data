%20160525 Xiaopang
% 比较 OA(72,4^1 3^24 2^20,2) 下，GMA-based designs 和 MWA-based designs 差异 
clear;clc;
epsilon = 1e-10; T0 = 1e-2;  T1 = 1e-6; 
InIter = 1e4; OutIter = 100*InIter; Reps = 30;
%InIter = 1e3; OutIter = 10*InIter; Reps = 3;
DiscMeasure = 'MD';

OA = importdata('OA from web/MA.72.4.1.3.24.2.20.txt');
[N,ncols] = size(OA);
q = zeros(ncols,1);
for i = 1:ncols
    q(i) = length(unique(OA(:,i)));
end
[q,id] = sort(q);
OA = OA(:,id);

Nsamp = 1000;
for n = 2:16
    q0 = [2*ones(n,1);3*ones(n,1)];
    aveDisc = zeros(Nsamp,1);
    A = zeros(Nsamp,2*n);
    WB = zeros(Nsamp,2*n);
    cols = zeros(Nsamp,2*n);
    for i = 1:Nsamp
        cols(i,:) = rand_subOA(OA,q0);
        D0 = OA(:,cols(i,:));
        [WB(i,:),aveDisc(i),A(i,:)] = Weighted_wordtype_pattern(D0,q0,DiscMeasure);
    end
    k = 10; %保留到小数点后第k位
    A = round(A*10^k)/10^k;
    WB = round(WB*10^k)/10^k;
    aveDisc = round(aveDisc*10^k)/10^k;
    [~,iwb] = sortrows(WB);
    [~,idisc] = sort(aveDisc);
    [~,ia] = sortrows(A);
    [~, Disc_vec0] = TA_MD_LHD(OA(:,cols(idisc(1),:)),q0,OutIter,InIter,T0,T1,Reps,0,1);
    [~, Disc_vec1] = TA_MD_LHD(OA(:,cols(ia(1),:)),q0,OutIter,InIter,T0,T1,Reps,0,1);
    fprintf('%d   %.6e  %.6e  %.6e  %.6e  %.6e  %.6e\n',n,mean(Disc_vec0),std(Disc_vec0),min(Disc_vec0),mean(Disc_vec1),std(Disc_vec1),min(Disc_vec1));
end
        
