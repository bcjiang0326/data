function test_SA_TA_SAIPM(q,n,isWD,choice,Reps,writedata)
% 20130205 测试算法 SA-REM,TA-REM,SA-IPM
% INPUT:
%       q: 各个因子水平数
%       n: 设计的行数
%       isWD: 0,CD; 1,WD
%       choice: 'SA', 'TA','SAIPM'
%       Reps: 重复次数 (default 1)
%       writedata: 0,不写文件(default); 1,写.
N = prod(q);

if nargin < 5
    Reps = 1;
    writedata = 0;
elseif nargin < 6
    writedata = 0;
end
if Reps > 1
    writedata = 0;
end

if  strcmp(choice,'SA')
    T0 = 1e-4;
    T1 = 1e-7;
    beta_in = 0.9;
    InIter = 10;
    OutIter = 500;
    t = cputime;
    if isWD
        [Disc,Disc_vec] = SA_REM_WD(n,q,OutIter,InIter,T0,T1,beta_in,Reps,writedata);
    else
        [Disc,Disc_vec] = SA_REM_CD(n,q,OutIter,InIter,T0,T1,beta_in,Reps,writedata);
    end
    t = cputime-t;
elseif strcmp(choice,'TA')
    Iter = 30000;
    InIter = 100;
    T0 = 1e-3;
    T1 = 1e-6;
    t = cputime;
    if isWD
        [Disc, Disc_vec] = TA_REM_WD(n,q,Iter,InIter,T0,T1,Reps,writedata);
    else
        [Disc, Disc_vec] = TA_REM_CD(n,q,Iter,InIter,T0,T1,Reps,writedata);
    end
    t = cputime-t;
elseif strcmp(choice,'SAIPM')
    if isWD
        % 实际上经过测试，WD情形似乎也应该采用CD的初始结论
        T0 = 1/q(1); % WD 情形
        beta_in = 0.99;
    else
        T0 = 1/(0.5*q(1)); %CD 情形
        beta_in = 0.99;
    end    
    InIter = 10;
    if N < 500
        OutIter = 10;
        beta_out = 0.9;
    else
        OutIter = 2;
        beta_out = 0.8;
    end
    t = cputime;
    [Disc,Disc_vec] = SA_IPM(n,q,OutIter,InIter,T0,beta_in,beta_out,isWD,Reps);
    t = cputime-t;
else
    error('Wrong choice!\n'); 
end

%D = Id2Design(q,ID);
%cout_8f([WD2_value(D,q),CD2_value(D,q)]);
cout_8f([Disc,mean(Disc_vec),std(Disc_vec),t]);

if writedata
    T_Disc_Disc0 = importdata('T_Disc_Disc0.txt');
    X1 = (1:length(T_Disc_Disc0))';
end

if strcmp(choice,'SA') && writedata
    plot(X1,T_Disc_Disc0(:,3),'LineWidth',0.01);%,X1,T_Disc_Disc0(:,4),'.','LineWidth',0.01);
elseif strcmp(choice,'TA') && writedata
    plot(X1,T_Disc_Disc0(:,3),'-','LineWidth',0.01); 
    %set(gca,'XTick',0:10000:100000)
    %set(gca,'XTickLabel',{'0','10,000','20,000','30,000','40,000',...
    %    '50,000','60,000','70,000','80,000','90,000','100,000'});
    %set(gca,'YtickLabel',{'0.1500','0.1550','0.1600','0.1650',...
    %    '0.1700','0.1750','0.1800','0.1850','0.1900'});
    %axis([0 100000 0.1500 0.1900])
    set(gca,'YGrid','on');
    %legend('D(96,3^54^5)',1);
    %ylabel('CD^2(D^c)');
elseif strcmp(choice,'SAIPM') && writedata
    plot(X1,T_Disc_Disc0(:,3),'-','LineWidth',0.01); 
end

end
