function [rate_wb,rate_a,Num] = comp_ranks(n1,n2,Measure)
%clear;
%n1 = 2; n2 = 2;
%Measure = 'MD';
D = importdata('OA from web/MA.36.3.12.2.11.txt');
D = D(:,[13:23,1:12]);

cols = enumerate_groups([11;12],[n1;n2]);
cols(:,n1+1:end) = cols(:,n1+1:end)+11;

q0 = [2*ones(n1,1);3*ones(n2,1)];

N = size(cols,1); n = n1+n2;

fid = fopen('out.txt','w');
for i = 1:N
    [WB,aveDisc,A,BB] = Weighted_wordtype_pattern(D(:,cols(i,:)),q0,Measure);
    fprintf(fid,'%.12e ',WB);
    fprintf(fid,'%.12e ',A);
    fprintf(fid,'%.12e ',aveDisc);
    fprintf(fid,'%.12e ',BB);
    fprintf(fid,'\n');
end
fclose(fid);

data = importdata('out.txt');
WB = data(:,1:n);
A = data(:,n+1:2*n);
aveDisc = data(:,2*n+1);
BB = data(:,2*n+2:end);

%设置精度
k = 10; %保留到小数点后第k位
A = round(A*10^k)/10^k;
BB = round(BB*10^k)/10^k;
WB = round(WB*10^k)/10^k;
aveDisc = round(aveDisc*10^k)/10^k;

%去掉字型型相同的行
[BB,ib] = unique(BB,'rows');
cols = cols(ib,:);
A = A(ib,:);
WB = WB(ib,:);
aveDisc = aveDisc(ib,:);

%按照平均均匀性排序
[aveDisc,id] = sort(aveDisc);
cols = cols(id,:);
A = A(id,:);
WB = WB(id,:);
BB = BB(id,:);
[~,ia] = sortrows(A);
[~,ia] = sort(ia);
[~,iwb] = sortrows(WB);
[~,iwb] = sort(iwb);

%计算p(r_A,r_disc),p(r_wb,r_disc)
Num = size(cols,1);
rate_a = 0;
rate_wb = 0;
for i = 1:Num-1
    for j = i+1:Num
        if ia(j)>ia(i)
            rate_a  = rate_a+1;
        end
        if iwb(j)>iwb(i)
            rate_wb = rate_wb+1;
        end
    end
end
rate_a = rate_a*2/(Num*(Num-1));
rate_wb = rate_wb*2/(Num*(Num-1));


