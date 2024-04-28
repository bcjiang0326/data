function [Disc0, Disc_vec, D0] = ...
    TA_adjusted_local_Disc(D0,Jset,OutIter,InIter,aa,cc,flag,RR,Reps,outfile)
%20190814 by Bochuan Jiang: adjusted TA algorithm
%参考Fang, Ke and Elsawah 2017 《Construction of uniform designs via
%    an adjusted threshold accepting algorithm》
%注意：本函数可适用于非对称设计
%INPUT:
%    D0: N-by-n matrix, the historical optimal design
%    Jset: k-by-1 vector, k <= n, 元素升序排列，标记D0中可改变的列
%    OutIter: 文章算法中的I，外层迭代次数，参考Example1，建议20
%    InIter: 文章算法中的J，内层迭代次数，参考Example1，建议5000
%    aa: 文章算法中的alpha, a scale or Reps-by-1 vector, 分别对应 section2.2
%        的 strategy 1 和 strategy 2, 参考Example1，建议0.15
%    cc: 文章算法中的c，参考Example2，建议0.03
%    flag:'CD','WD','MD'
%    RR: RR=range(disc)，用于计算 T1=aa*RR,
%    Reps: 重复次数 (default,1)
%    outfile: 文件名，将最后一次重复的结果写入文件 (default,'');
%OUTPUT:
%    Disc0: 得到的最优设计的Disc值
%    Disc_vec: Reps-by-1 vector, 记录每次重复的最优解
%    D0: 得到的最优矩阵

epsilon = 1e-9;
[N,n] = size(D0);
q = zeros(n,1);
for i = 1:n
    q(i) = length(unique(D0(:,i)));
end
if ~issorted(Jset) || any(Jset~=floor(Jset)) || max(Jset) > n || min(Jset) < 1 ||...
        length(unique(Jset))~=length(Jset) || size(Jset,2) ~= 1 || length(size(Jset))>2
    error('Wrong Jset!\n');
end
if OutIter < 1 || InIter < 1
    error('Wrong OutIter, InIter!\n');
end
if ~(0 < cc && cc < 1)
    error('Wrong cc!\n');
end
if nargin < 7 || isempty(flag)
    flag = 'CD';
end
if ~strcmp(flag,'CD') && ~strcmp(flag,'WD') && ~strcmp(flag,'MD')
    error('Wrong flag!\n');
end
if nargin < 8 || isempty(RR)
    RR = Disc_Range(D0,Jset,InIter,flag);
end
if RR < 0
    error('Wrong RR!\n');
end
if nargin < 9 || isempty(Reps) || Reps < 1
    Reps = 1;
end
if length(size(aa))>2 || size(aa,2)~=1 || any(aa<=0) || any(aa>1) ||...
        ( size(aa,1)>1 && size(aa,1)~=Reps )
    error('Wrong aa!\n');
end
if size(aa,1)==1
    aa = aa*ones(Reps,1);
end
if nargin < 10 || isempty(outfile)
    iswrite = 0;
else
    iswrite = 1; % 设置输出
    fid = fopen(outfile,'w');
end

%首先计算阈值序列
T = zeros(OutIter,Reps);
T(1,:) = RR*aa';
for rep = 1:Reps
    for t = 2:OutIter-1
        if T(t-1,rep) > cc*T(1,rep)
            T(t,rep) = (OutIter-t)/OutIter*T(t-1,rep);
        else
            T(t,rep) = (OutIter-t-1)/(OutIter-t)*T(t-1,rep);
        end
    end
end


%计算给定的 historical optimal design 的sigma0,sigma01,Disc0和z0
if strcmp(flag,'MD')
    [Disc0,sigma0,sigma01] = MD2_value(D0,q);
    z0 = (D0+0.5)./(ones(N,1)*q')-0.5;
elseif strcmp(flag,'CD')
    [Disc0,sigma0,sigma01] = CD2_value(D0,q);
    z0 = 0.5*( (D0+0.5)./(ones(N,1)*q')-0.5 );
elseif strcmp(flag,'WD')
    [Disc0,sigma0] = WD2_value(D0,q);
end

%记录历次重复最后得到的最优解
Disc_vec = zeros(Reps,1);

%进入算法
D = D0; %current design
for rep = 1:Reps
    % 首先对Jset标记的列进行随机化，作为第rep重复的初始矩阵
    for j = Jset'
        D(:,j) = D0(randperm(N,N),j);
    end
    % 生成初始的sigma,sigma1,Disc 和 z
    if strcmp(flag,'MD')
        [Disc,sigma,sigma1] = MD2_value(D,q);
        z = (D+0.5)./(ones(N,1)*q')-0.5;
    elseif strcmp(flag,'CD')
        [Disc,sigma,sigma1] = CD2_value(D,q);
        z = 0.5*( (D+0.5)./(ones(N,1)*q')-0.5 );
    elseif strcmp(flag,'WD')
        [Disc,sigma] = WD2_value(D,q);
    end
    
    %当次重复的初始解在全部重复中更优？是，更新D0,sigma0,sigma01
    Disc_vec(rep) = Disc;
    if Disc_vec(rep) < Disc0-epsilon
        D0 = D;
        Disc0 = Disc_vec(rep);
        sigma0 = sigma;
        if strcmp(flag,'CD') || strcmp(flag,'MD')
            sigma01 = sigma1;
            z0 = z;
        end
    end
    
    %进入外层循环
    for Oiter = 1:OutIter
        %Current Design差于D0？是，取D0作为Current Design
        if Disc > Disc0+epsilon
            D = D0;
            Disc = Disc0;
            sigma = sigma0;
            if strcmp(flag,'CD') || strcmp(flag,'MD')
                sigma1 = sigma01;
                z = z0;
            end
            %新的Current Design在本次重复中更优？是，更新本次重复最优解
            if Disc < Disc_vec(rep) - epsilon
                Disc_vec(rep) = Disc;
            end
        end
        
        %进入内层循环
        for Iiter = 1:InIter
            %产生要交换的列及列中两个元素
            k = Jset(randi(length(Jset)));
            i = randi(N);
            j = randi(N);
            while j==i || abs(D(i,k)-D(j,k))<epsilon
                j=randi(N);
            end
            
            %产生用于更新 delta 的 alpha1,alpha2,beta
            if strcmp(flag,'MD')
                alpha1=(5/3-0.25*abs(z(j,k))-0.25*z(j,k)^2)/(5/3-0.25*abs(z(i,k))-0.25*z(i,k)^2);
                alpha2=(1.875-0.5*abs(z(j,k)))/(1.875-0.5*abs(z(i,k)));
                beta = (1.875-0.25*abs(z(j,k))-0.25*abs(z(:,k))-0.75*abs(z(j,k)-z(:,k))+0.5*(z(j,k)-z(:,k)).^2)...
                    ./(1.875-0.25*abs(z(i,k))-0.25*abs(z(:,k))-0.75*abs(z(i,k)-z(:,k))+0.5*(z(i,k)-z(:,k)).^2);
            elseif strcmp(flag,'CD')
                alpha1=(1+abs(z(j,k))-2*z(j,k)^2)/(1+abs(z(i,k))-2*z(i,k)^2);
                alpha2=(1+2*abs(z(j,k)))/(1+2*abs(z(i,k)));
                beta = (1+abs(z(j,k))+abs(z(:,k))-abs(z(j,k)-z(:,k)))...
                    ./(1+abs(z(i,k))+abs(z(:,k))-abs(z(i,k)-z(:,k)));
            elseif strcmp(flag,'WD')
                tempj = abs(D(j,k)-D(:,k))/q(k);
                tempi = abs(D(i,k)-D(:,k))/q(k);
                beta = ( 1.5-tempj.*(1-tempj) )./( 1.5-tempi.*(1-tempi) );
            end
            
            %计算delta
            delta = 0;
            for t = 1:N
                if t ~= i && t ~= j
                    delta = delta ...
                        +(beta(t)-1)*(sigma(i,t)-sigma(j,t)/beta(t));
                end
            end
            delta = 2*delta;
            if strcmp(flag,'CD') || strcmp(flag,'MD')
                delta = delta +(alpha1-1)*(sigma1(j)/alpha1-sigma1(i))...
                    +(alpha2-1)*(sigma(i,i)-sigma(j,j)/alpha2);
            end
            
            %判断是否更新
            if delta < T(Oiter,rep)-epsilon
                % 更新 D,Disc
                temp=D(i,k);D(i,k)=D(j,k);D(j,k)=temp;
                Disc = Disc + delta;
                % 更新 sigma 非对角部分
                for t=1:N
                    if t~=i && t~=j
                        sigma(i,t)=sigma(i,t)*beta(t);
                        sigma(j,t)=sigma(j,t)/beta(t);
                        sigma(t,i)=sigma(i,t);
                        sigma(t,j)=sigma(j,t);
                    end
                end
                if strcmp(flag,'CD') || strcmp(flag,'MD')
                    % 更新 z,sigma1,diag(sigma)
                    temp=z(i,k);z(i,k)=z(j,k);z(j,k)=temp;
                    sigma1(i)=sigma1(i)*alpha1;
                    sigma1(j)=sigma1(j)/alpha1;
                    sigma(i,i)=sigma(i,i)*alpha2;
                    sigma(j,j)=sigma(j,j)/alpha2;
                end
                
                % 当 delta<0 说明开始下降.
                if delta < -epsilon
                    % 若遇到比本次更优解，则记录之
                    if Disc < Disc_vec(rep) - epsilon
                        Disc_vec(rep) = Disc;
                        % 若遇到全局更优解，则记录之
                        if Disc < Disc0-epsilon
                            D0 = D;
                            Disc0 = Disc;
                            sigma0 = sigma;
                            if strcmp(flag,'CD') || strcmp(flag,'MD')
                                sigma01 = sigma1;
                                z0 = z;
                            end
                        end
                    end
                end
            end
            if iswrite && rep == Reps
                fprintf(fid,'%.8f\n',Disc);
            end
        end
    end
end

if iswrite
    fclose(fid);
end
end


%{
%测试代码
N = 27; n = 13; s = 3;
q = s*ones(n,1);
D0 = rand_U_Type_design(q,N);
Jset = (1:n)';

InIter = 5000;
OutIter = 20;
aa = [0.15,0.016,0.01,0.002,0.0005]';
cc = 0.03;
flag = 'WD';
Reps = 1;
outfile = 'out.txt';

figure
for i = 1:length(aa)
    [Disc1, Disc_vec, D0] = TA_adjusted_local_Disc(D0,Jset,OutIter,InIter,aa(i),cc,flag,'',Reps,outfile);
    if strcmp(flag,'MD')
        Disc = MD2_value(D0);
    elseif strcmp(flag,'CD')
        Disc = CD2_value(D0);
    elseif strcmp(flag,'WD')
        Disc = WD2_value(D0);
    end
    fprintf('%.6f  %.6f  %.6f  %d\n', Disc1, Disc, mean(Disc_vec),isOA(D0,2));
    aveDisc_seq = importdata('out.txt');    
    plot(1:length(aveDisc_seq),aveDisc_seq);
    hold on
end
hold off
%}
