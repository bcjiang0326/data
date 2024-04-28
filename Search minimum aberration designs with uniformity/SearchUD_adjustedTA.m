flag = 'WD';
Reps = 30;
InIter = 10000; %�ڲ�����ͳͳ����Ϊ1��

%Ԫ���û��ĵ�������
OutIter1 = 100;
AllIter1 = OutIter1*InIter;

%ˮƽ�û��ĵ�������
OutIter2 = 10;
AllIter2 = OutIter2*InIter;

epsilon = 1e-8;
oa = importdata('OA from web/oa.32.9.4.2.a.txt');
[N,n_oa] = size(oa);
s = length(unique(oa(:,1)));

%�㷨����
aa = 0.15;
cc = 0.03;
T0 = 1e-2; 
T1 = 1e-6;

delta = 1;
for n = 53:53
    q = s*ones(n,1);
    
    Jset2 = (n_oa+1:n)';
    D0 = [ oa, rand_U_Type_design(s*ones(n-n_oa,1),N) ];
    %��ʼ��ʱ
    t2 = cputime;
    %����2����1. Ԫ���û�����
    RR = Disc_Range(D0,Jset2,InIter,flag);
    [Disc21,~,D1] = TA_adjusted_local_Disc(D0,Jset2,OutIter1,InIter,aa,cc,flag,RR,Reps);
    %����2����2. ˮƽ�û�����
    [Disc2, Disc2vec, D2] = TA_uniform_FF_design_WD(D1,AllIter2,InIter,T0,T1,Reps);
    
    t2 = cputime-t2;
    A2 = GWP(D2);
    ProjDisc2 = Proj_Disc2(2, D2, flag);%��άͶӰ������
    
    %д�ļ�
    file = strcat('My UD/',flag,'/level',int2str(s),'/UD.',int2str(N),'.',int2str(s),'.',int2str(n),'.txt');
    Disc_ex = inf;
    if exist(file,'file')
        D_ex = importdata(file);
        Disc_ex = Disc2_value(D_ex,'',flag);
        A_ex = GWP(D_ex);
        ProjDisc_ex = Proj_Disc2(2, D_ex, flag);
    end
    if Disc2 < Disc_ex * (1-epsilon)
        fid2 = fopen(file,'w');
        for row = 1:size(D2,1)
            fprintf(fid2,'%d ',D2(row,1:end-1));
            fprintf(fid2,'%d\n',D2(row,end));
        end
        fclose(fid2);
        Disc_ex = Disc2;
        A_ex = A2;
        ProjDisc_ex = ProjDisc2;
    end
    
    y = [n,t2,Disc_ex,ProjDisc_ex,A_ex(2),Disc2,ProjDisc2,A2(2)];
    fprintf('%d  %.1f | %.6f   %.6f   %.1f | %.6f   %.6f   %.1f   \n',y);
end
