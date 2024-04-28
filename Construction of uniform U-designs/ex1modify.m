%20151022 by Xiaopang �����������°��� MD ׼������
clear;
load ex1.mat
q0 = [2;2;2;3;3;3];
k = 13; %������С������kλ
A = round(A*10^k)/10^k;
BB = round(BB*10^k)/10^k;

[BB,ib] = unique(BB,'rows');
cols = cols(ib,:);
A = A(ib,:);

aveMD = zeros(size(BB,1),1); % ����� U-design ��ƽ�� MD ֵ
WB_MD = zeros(size(BB,1),length(q0)); % weighted wordtype pattern

%aveCD = zeros(size(BB,1),1); % ����� U_design ��ƽ�� CD ֵ
%WB_CD = zeros(size(BB,1),length(q0)); % weighted wordtype pattern

%aveWD = zeros(size(BB,1),1); % ����� U_design ��ƽ�� WD ֵ
%WB_WD = zeros(size(BB,1),length(q0)); % weighted wordtype pattern

for i = 1:size(cols,1)
    [WB_MD(i,:),aveMD(i)] = Weighted_wordtype_pattern(D(:,cols(i,:)),q0,'WD');
end
WB_MD = WB_MD(:,[3:end]); % ����ǿ���� 2, ��WB_MD(1)=WB_MD(2)=0;
WB_MD = round(WB_MD*10^k)/10^k; % ��������ΪС�����12λ
aveMD = round(aveMD*10^k)/10^k; % ��������ΪС�����12λ

%{
for i = 1:size(cols,1)
    [WB_CD(i,:),aveCD(i)] = Weighted_wordtype_pattern(D(:,cols(i,:)),q0,'CD');
end
WB_CD = WB_CD(:,[3:end]); % ����ǿ���� 2, ��WB_CD(1)=WB_CD(2)=0;
WB_CD = round(WB_CD*10^k)/10^k; % ��������ΪС�����12λ
aveCD = round(aveCD*10^k)/10^k; % ��������ΪС�����12λ

for i = 1:size(cols,1)
    [WB_WD(i,:),aveWD(i)] = Weighted_wordtype_pattern(D(:,cols(i,:)),q0,'WD');
end
WB_WD = WB_WD(:,[3:end]); % ����ǿ���� 2, ��WB_WD(1)=WB_WD(2)=0;
WB_WD = round(WB_WD*10^k)/10^k; % ��������ΪС�����12λ
aveWD = round(aveWD*10^k)/10^k; % ��������ΪС�����12λ
%}

%����������ƵĽ������ MD ֵ����
[aveMD,imd] = sort(aveMD);
cols = cols(imd,:);
BB = BB(imd,:);
A = A(imd,:);
WB_MD = WB_MD(imd,:);

%aveCD = aveCD(imd);
%WB_CD = WB_CD(imd,:);

%aveWD = aveWD(imd);
%WB_WD = WB_WD(imd,:);


[~,ia] = sortrows(A);
[~,ia] = sort(ia);
[~,iwbmd] = sortrows(WB_MD);
[~,iwbmd] = sort(iwbmd);

%[~,icd] = sort(aveCD);
%[~,iwbcd] = sortrows(WB_CD);

%[~,iwd] = sort(aveWD);
%[~,iwbwd] = sortrows(WB_WD);

%
Num = 20;
for i = 1:Num
   x = [cols(i,:),aveMD(i),i,iwbmd(i),ia(i)]; 
   fprintf('%d %d %d %d %d %d & %.6f & %d & %d & %d \\\\ \n',x);
   %x = [cols(i,:),i,iwbmd(i),icd(i),iwbcd(i),iwd(i),iwbwd(i),ia(i)];
   %fprintf('%d %d %d %d %d %d & %d & %d & %d & %d & %d & %d & %d \\\\ \n',x);
end
%}


%{
%��ͼ
Disc_fun = @MD2_value;
k1 = 1; k2 = 674; 
D1 = D(:,cols(k1,:));
D2 = D(:,cols(k2,:));

[N,n] = size(D1);
q_lh = N*ones(n,1);
K = 1e6; 
vec1 = zeros(K,1); 
vec2 = zeros(K,1);
for i = 1:K
    L1 = rand_oa_lhd(D1,q0,1);
    vec1(i) = Disc_fun(L1,q_lh);
    L2 = rand_oa_lhd(D2,q0,1);
    vec2(i) = Disc_fun(L2,q_lh);
end

ymin=min([vec1;vec2]);
ymax=max([vec1;vec2]);
x=linspace(ymin,ymax,50); %�������С����ֳ�20���ȷֵ�(19�ȷ�),Ȼ��ֱ�����������ĸ��� 
yy1=hist(vec1,x);  %�����������ĸ���
yy2=hist(vec2,x);  %�����������ĸ���

figure;
plot(x,yy1,'k-o',x,yy2,'k-*'); hold on;
xlim([ymin,ymax]);
bar(x,[yy1;yy2]','grouped','k'); 
legend('D_1-based','D_2-based');
hold off;
%}


