function [Disc0, Disc_vec, L0] = TA_MD_LHD_1(groups,D,q,OutIter,InIter,T0,T1,Reps,writedata)
% 20150126 
% 基于 TA 算法求解均匀 LHD
% 求解均匀设计，迭代采用 columnwise-pairwise 策略
% 背景：Tang, Xu and Lin (2012) AOS 提出了uniform fractional factorial design (UFFD). 基于此
%       Tang and Xu (2012) JSPI 提出了一种构造 uniform FFD 的启发算法。 本函数基于其构造的
%       UFFD 构造 OA-based LHD. 构造过程采用 TA 算法。
% 准则：Zhou, Fang and Ning (2013) Mixture discrepancy
% 使用范围：根据给定的D，计算 D-based LHD，这里仅要求 D 是 balanced。 
% Neighborhood: 记 D, ( N-by-s design) 为 UFFD，D 的第一列为 q1 水平。令 k = n/q1，采用如下变换
%       0 ---> {0,...,k-1}， 1 ---> {k,...,2k-1}, 依次类推。将所有列均如此变换得到一个 LHD 记为 L。
%       此时 L 的 Neighborhood 定义为 置换 L 列中由 D 的同一元素对应的行得到的所有 LHD 的集合
%       .
% INPUT:
%       D: N-by-n balanced design
%       q: n-by-1 vector, D 各因子的水平数
%       OutIter: 外部迭代次数
%       InIter: 内部迭代次数
%       T0: 初始阈值，位于(0,1)之间，建议1e-2
%       T1: 最终阈值，位于(0,1)之间，建议1e-6
%       Reps: 重复次数 (default,1)
%       writedata: 是否将过程数据写入文件。 1,写入; 0(default) 不写
% OUTPUT:
%       Disc0: 得到最优矩阵的 Squared Discrepancy
%       Disc_vec: Reps-by-1 vector, 记录每次重复的最优解
%       L0: 得到的最优矩阵

[N,n] = size(D);
if size(q,2)~=1 || n~=size(q,1)
    error(' Input variable q must be a s-by-1 vector!\n');
end
if any(N./q~=floor(N./q))
    error('N and q does not match each other !\n');
end
if nargin < 7
    Reps = 1;
    writedata = 0;
elseif nargin < 8
    writedata = 0;
end
if Reps > 1
    writedata = 0;
end

% 设置文件输出
if writedata
    fid = fopen('T_Disc_Disc0.txt','w');
end

% 设置归零参数
epsilon = 1e-12;

% 首先计算阈值变化的比率,即每 InIter 次迭代，
% 更新一次阈值 T = T*ratio     
ratio = (T1/T0)^(1/(ceil(OutIter/InIter)-1));

% 计算初始的 UFFD-based LHD。首先，依次记录各列的 group 划分情况
% groups 第一列前 n/q1 个元素对应 L 第一列的第一个 group，依次类推。
KK = N./q; % 记录各列中，每个group 元素个数
L = zeros(size(D));

Disc0 = inf; 
Disc_vec = zeros(Reps,1);

for rep = 1:Reps
    % 随机产生一个 L
    for k = 1:n
        for i = 1:q(k)
            v = randperm(KK(k));
            L( groups((i-1)*KK(k)+v,k), k ) = (0:KK(k)-1)'+(i-1)*KK(k);
        end
    end
    
    % 首先计算初始的 z,sigma1,sigma
    z = (L+0.5)/N-0.5;
    sigma1 = zeros(N,1);
    sigma = zeros(N,N);
    for i = 1:N-1
        for j = i+1:N
            sigma(i,j) = 1/N^2;
            for k = 1:n
                sigma(i,j)=sigma(i,j)*(1.875-0.25*abs(z(i,k))-0.25*abs(z(j,k))-0.75*abs(z(i,k)-z(j,k))+0.5*(z(i,k)-z(j,k))^2);
            end
            sigma(j,i)=sigma(i,j);
        end 
    end
    for i = 1:N
        sigma1(i) = 2/N;
        for k = 1:n
            sigma1(i) = sigma1(i)*(5/3-0.25*abs(z(i,k))-0.25*z(i,k)^2);
        end
        sigma(i,i) = prod(1.875-0.5*abs(z(i,:)))/N^2;
    end
    
    % 初始化Disc
    Disc=0;
    for i = 1:N-1
        for j = i+1:N
            Disc = Disc + sigma(i,j);
        end
    end
    Disc=(19/12)^n+2*Disc+sum(diag(sigma))-sum(sigma1);

    % 判断当次重复的初始解是否更优
    Disc_vec(rep) = Disc; 
    if Disc_vec(rep) < Disc0-epsilon
        Disc0 = Disc_vec(rep);
        if nargout > 2
            L0=L;
        end
    end
    
    %进入算法
    T = T0;
    for Oiter = 1:OutIter
        if mod(Oiter,InIter)==0
            T = T*ratio;
        end
        %产生要交换的列及列中两个元素    
        k=randi(n); % 随机选第 k 列
        i0 = randi(q(k)); %随机选该列中的第 i0 个group
        i1=randi( KK(k) ); % 随机选该 group 中第 i1 和 i2 个元素
        i2=randi( KK(k) );
        while i2==i1    
            i2=randi( KK(k) );
        end
        i = groups( (i0-1)*KK(k)+i1, k);
        j = groups( (i0-1)*KK(k)+i2, k);
        if D(i,k)~=D(j,k)
            error('shit!\n');
        end
    
        %产生用于计算 delta 的 alpha1,alpha2,beta
        alpha1=(5/3-0.25*abs(z(j,k))-0.25*z(j,k)^2)/(5/3-0.25*abs(z(i,k))-0.25*z(i,k)^2);
        alpha2=(1.875-0.5*abs(z(j,k)))/(1.875-0.5*abs(z(i,k)));
        beta = (1.875-0.25*abs(z(j,k))-0.25*abs(z(:,k))-0.75*abs(z(j,k)-z(:,k))+0.5*(z(j,k)-z(:,k)).^2)...
            ./(1.875-0.25*abs(z(i,k))-0.25*abs(z(:,k))-0.75*abs(z(i,k)-z(:,k))+0.5*(z(i,k)-z(:,k)).^2);
    
        %计算delta          
        delta = 0;
        for t = 1:N
            if t ~= i && t ~= j
                delta = delta ...
                    +(beta(t)-1)*(sigma(i,t)-sigma(j,t)/beta(t));
            end
        end
        delta = 2*delta+(alpha1-1)*(sigma1(j)/alpha1-sigma1(i))...
            +(alpha2-1)*(sigma(i,i)-sigma(j,j)/alpha2);
        
        %判断是否更新
        if delta <  Disc*T-epsilon
        %if  T < 1e-5 || delta <  Disc*T-epsilon
            % 更新 L,Disc
            temp=L(i,k);L(i,k)=L(j,k);L(j,k)=temp;
            Disc = Disc + delta;
            % 更新 z,sigma1,diag(sigma),sigma
            temp=z(i,k);z(i,k)=z(j,k);z(j,k)=temp;            
            sigma1(i)=sigma1(i)*alpha1;
            sigma1(j)=sigma1(j)/alpha1;
            sigma(i,i)=sigma(i,i)*alpha2;
            sigma(j,j)=sigma(j,j)/alpha2;
            for t=1:N
                if t~=i && t~=j
                    sigma(i,t)=sigma(i,t)*beta(t);
                    sigma(j,t)=sigma(j,t)/beta(t);
                    sigma(t,i)=sigma(i,t);
                    sigma(t,j)=sigma(j,t);
                end
            end
            % 当 delta<0 说明开始下降.
            if delta < -epsilon
                % 若遇到比本次更优解，则记录之
                if Disc < Disc_vec(rep) - epsilon
                    Disc_vec(rep) = Disc;
                    % 若遇到全局更优解，则记录之
                    if Disc_vec(rep) < Disc0-epsilon
                        Disc0 = Disc_vec(rep);
                        if nargout > 2
                            L0 = L;
                        end
                    end
                end
            end
        end
        if writedata
            fprintf(fid,'%.8f %.8f %.8f %.8f %.8f\n',delta/Disc,T,Disc,Disc0,Oiter);
        end
    end
end

if writedata
    fclose(fid);
end
end

%{
% 测试函数
N = 27; n = 5; q = 3*ones(n,1);
X = combins(3,3); % 自由列
C = [1 1 1; 1 2 0; 1 1 2; 1 0 1; 0 1 2; 1 2 2; 1 1 0]';
B = { [0,1], [0,1,2], [1,1,0,2], [1,1,0,2,1], [1,1,0,2,1,2], [1,1,0,2,1,2,2] };
D0 = mod( [X,X*C(:,1:n-3)+ones(N,1)*B{n-4}], 3);
InIter = 100;
OutIter = 100*InIter;
Reps = 3;
[Disc0, Disc_vec, L0] = TA_MD_LHD(D0,q,OutIter,InIter,1e-3,1e-6,1,0);
fprintf('%.8f\n',Disc0);
fprintf('%.8f\n',MD2_value(L0,ones(n,1)*N));
%}



