%生成例1中的表2

OA = importdata('OA from web/MA.36.3.12.2.11.txt');
[N,ncols] = size(OA);
q = zeros(ncols,1);
for i = 1:ncols
    q(i) = length(unique(OA(:,i)));
end
[q,id] = sort(q);
OA = OA(:,id);

n10 = 3;
n20 = 3;
n = n10+n20;
q_oa = [2*ones(n10,1); 3*ones(n20,1)];
q_lh = N*ones(n,1);

filename = 'MA LHD 20151102 MD oa36 comp order/';
filename = strcat(filename,'results',int2str(n10),int2str(n20));
    
results = importdata(filename);
cols = results(:,1:n);
Num = size(cols,1);

aveMD = zeros(Num,1);
aveCD = zeros(Num,1);
aveWD = zeros(Num,1);
WBmd = zeros(Num,n);
WBcd = zeros(Num,n);
WBwd = zeros(Num,n);
A = zeros(Num,n);

for i = 1:Num
    if mod(i,100)==0
        fprintf('%d\n',i);
    end
    [WBmd(i,:),aveMD(i),A(i,:)] = Weighted_wordtype_pattern(OA(:,cols(i,:)),q_oa,'MD');
    [WBcd(i,:),aveCD(i)] = Weighted_wordtype_pattern(OA(:,cols(i,:)),q_oa,'CD');
    [WBwd(i,:),aveWD(i)] = Weighted_wordtype_pattern(OA(:,cols(i,:)),q_oa,'WD');
end

k = 12;
aveMD = round(aveMD*10^k)/10^k; % 保留精度为小数点后k位
aveCD = round(aveCD*10^k)/10^k; % 保留精度为小数点后k位
aveWD = round(aveWD*10^k)/10^k; % 保留精度为小数点后k位
WBmd = round(WBmd*10^k)/10^k; % 保留精度为小数点后k位
WBcd = round(WBcd*10^k)/10^k; % 保留精度为小数点后k位
WBwd = round(WBwd*10^k)/10^k; % 保留精度为小数点后k位
A = round(A*10^k)/10^k; % 保留精度为小数点后k位

[aveMD,a] = sort(aveMD);
cols = cols(a,:);
aveCD = aveCD(a);
aveWD = aveWD(a);
WBmd = WBmd(a,:);
WBcd = WBcd(a,:);
WBwd = WBwd(a,:);
A = A(a,:);

[~,icd] = sort(aveCD);
[~,icd] = sort(icd);
[~,iwd] = sort(aveWD);
[~,iwd] = sort(iwd);
[~,iwb_md] = sortrows(WBmd);
[~,iwb_md] = sort(iwb_md);
[~,iwb_cd] = sortrows(WBcd);
[~,iwb_cd] = sort(iwb_cd);
[~,iwb_wd] = sortrows(WBwd);
[~,iwb_wd] = sort(iwb_wd);
[~,ia] = sortrows(A);
[~,ia] = sort(ia);
for i = 1:20
    fprintf('%d ',cols(i,:));
    fprintf('& %.6f & ',aveMD(i));
    y = [i,icd(i),iwd(i),iwb_md(i),iwb_cd(i),iwb_wd(i),ia(i)];
    fprintf('%d & %d & %d & %d & %d & %d & %d \\\\\n',y);
end


