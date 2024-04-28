% 用于测试各方法的最优参数
% SA_IPM: 选择 T0=1/(k*levels) 中的 k,及 beta_in
% SA_REM: 选择 T0, T1, beta_in
% Case 1
levels = 4;s = 5;n = 200; 
q = levels*ones(s,1);
N = levels^s;



choice = 0; %0,SA; 1,TA; 2,SAIPM
isWD = 0;
writedata = 0;
Reps = 100;
if  choice == 0
    %k1_zone = [2,3,4];
    %k2_zone = [4,5,6,7];
    k1_zone = [4];
    k2_zone = [6,7,8];
    Tf = 0.8;
    InIter = 10;
    OutIter = 500;
    for k1 = k1_zone
        for k2 = k2_zone
            T0 = 10^(-k1);
            T1 = 10^(-k2);
            cout_d([k1,k2]);
            t = cputime;
            if isWD
                [Disc,Disc_vec] = SA_REM_WD(n,q,OutIter,InIter,T0,T1,Tf,Reps,writedata);
            else
                [Disc,Disc_vec] = SA_REM_CD(n,q,OutIter,InIter,T0,T1,Tf,Reps,writedata);
            end
            cout_8f([Disc,mean(Disc_vec),cputime-t]); 
        end
    end
elseif choice == 1
    Iter = 30000;
    InIter = 100;
    T0 = 1e-3;
    T1 = 1e-6;
    t = cputime;
    if isWD
        [Disc, ID, D] = TA_REM_WD(n,q,Iter,InIter,T0,T1,writedata);
    else
        [Disc, ID, D] = TA_REM_WD(n,q,Iter,InIter,T0,T1,writedata);
    end
    cout_8f([Disc,cputime-t]);
elseif choice == 2
    k = 0.5;
    beta_in= 0.9;
    beta_out_zone = [0.95,0.9,0.85,0.8];
    %for k = k_zone
        for beta_out = beta_out_zone
            T0 = 1/(levels*k);
            InIter = 10;
            Reps = 50;
            if N < 500
                OutIter = 10;
                beta_out = 0.9;
            else
                OutIter = 2;
                beta_out = 0.8;
            end
            cout_2f(beta_out);
            t1 = cputime;
            [Disc,Disc_vec] = SA_IPM(n,q,OutIter,InIter,T0,beta_in,beta_out,isWD,Reps);
            t1 = cputime-t1;
            cout_8f([Disc,mean(Disc_vec),std(Disc_vec),t1]);
        end
    %end
end

if writedata
    T_Disc_Disc0 = importdata('T_Disc_Disc0.txt');
    X1 = (1:length(T_Disc_Disc0))';
end

if ~choice && writedata
    plot(X1,T_Disc_Disc0(:,3),'LineWidth',0.01);%,X1,T_Disc_Disc0(:,4),'.','LineWidth',0.01);
elseif choice==1 && writedata
    plot(X1,T_Disc_Disc0(:,3),'-','LineWidth',0.01); 
    %set(gca,'XTick',0:10000:100000)
    %set(gca,'XTickLabel',{'0','10,000','20,000','30,000','40,000',...
    %    '50,000','60,000','70,000','80,000','90,000','100,000'});
    %set(gca,'YtickLabel',{'0.1500','0.1550','0.1600','0.1650',...
    %    '0.1700','0.1750','0.1800','0.1850','0.1900'});
    %axis([0 100000 0.1500 0.1900])
    set(gca,'YGrid','on');
    legend('D(96,3^54^5)',1);
    ylabel('CD^2(D^c)');
elseif choice==2 && writedata
    plot(X1,T_Disc_Disc0(:,3),'-','LineWidth',0.01); 
end
