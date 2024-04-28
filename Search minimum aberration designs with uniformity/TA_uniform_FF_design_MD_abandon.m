function [MD0, MDvec, D0] = TA_uniform_FF_design_MD_abandon(D,q,OutIter,InIter,T0,T1,Reps,writedata)
% 20181219
% 采用 MD 准则
% 基于 TA 算法求解 uniform fractional factorial design (UFFD)
% 求解均匀设计，迭代采用 columnwise-pairwise 策略, 每次迭代中选择一列，选择两个水平，i,j
%       将水平为 i 的元素变为 j, 将水平为 j 的变为 i.
% 背景：Tang, Xu and Lin (2012) AOS 提出了uniform fractional factorial design (UFFD). 基于此
%       Tang and Xu (2012) JSPI 提出了一种构造 uniform FFD 的启发算法。 
% 使用范围：根据给定的D，计算所有水平置换下的均匀设计
% Neighborhood: 记 D 邻域为：D的任意一列两个水平置换得到的设计的集合
%       .
% INPUT:
%       D: N-by-n balanced design
%       q: s-by-1 vector, 设计各因子的水平数
%       OutIter: 外部迭代次数
%       InIter: 内部迭代次数
%       T0: 初始阈值，位于(0,1)之间，建议1e-2
%       T1: 最终阈值，位于(0,1)之间，建议1e-6
%       Reps: 重复次数 (default,1)
%       writedata: 是否将过程数据写入文件。 1,写入; 0(default) 不写
% OUTPUT:
%       MD0: 得到最优矩阵的 Squared Discrepancy
%       MD_vec: Reps-by-1 vector, 记录每次重复的最优解
%       D0: 得到的最优矩阵

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
    fid = fopen('T_Disc_Disc.txt','w');
end

% 设置归零参数
epsilon = 1e-10;

% 首先计算阈值变化的比率,即每 InIter 次迭代，
% 更新一次阈值 T = T*ratio     
ratio = (T1/T0)^(1/(ceil(OutIter/InIter)-1));

KK = N./q; % 记录各列中，每个水平重复次数

MD0 = inf; 
MDvec = zeros(Reps,1);

for rep = 1:Reps
    % 首先生成一个随机无重复的矩阵，生成初始的 z,sigma1,sigma2,sigma
    sigma1 = zeros(N,1);
    sigma2 = zeros(N,1);
    sigma = zeros(N,N);
    z = 0.25*( (D+0.5)./(ones(N,1)*q')-0.5);
    for i = 1:N-1
        for j = i+1:N
            sigma(i,j) = 1/N^2*prod(1.875-abs(z(i,:))-abs(z(j,:))-3*abs(z(i,:)-z(j,:))+8*(z(i,:)-z(j,:)).^2);
            sigma(j,i)=sigma(i,j);
        end 
    end
    for i = 1:N
        sigma1(i) = 2/N*prod(5/3-abs(z(i,:))-4*z(i,:).^2);
        sigma2(i) = 1/N^2*prod(1.875-2*abs(z(i,:)));
    end
    
    % 初始化Disc
    Disc=0;
    for i = 1:N-1
        Disc = Disc + sum(sigma(i,i+1:N));
    end
    Disc=(19/12)^n+2*Disc+sum(sigma2)-sum(sigma1);

    % 判断当次重复的初始解是否更优
    MDvec(rep) = Disc; 
    if MDvec(rep) < MD0-epsilon
        MD0 = MDvec(rep);
        D0=D;
    end
    
    %进入算法
    T = T0;
    for Oiter = 1:OutIter
        if mod(Oiter,InIter)==0
            T = T*ratio;
        end
        %产生要交换的列
        k=randi(n);
        %产生要交换的水平集
        i=randi(N);
        j=randi(N);
        while j==i || D(i,k)==D(j,k)    
            j=randi(N);
        end
        %置换2水平因子的0和1或者3水平因子的0和2，均不改变CD值.
        if q(k)==2 || ( q(k)==3 && D(i,k)+D(j,k) == 2 )
            if writedata
                fprintf(fid,'%.8f %.8f %.8f %.8f %.8f\n',0,T,Disc,MD0,Oiter);
            end
            continue;
        end
        Iset = D(:,k)==D(i,k); 
        Jset = D(:,k)==D(j,k);
        Tset = ~( Iset | Jset );
      
        %产生用于计算 delta 的 alpha1,alpha2,beta
        alpha1=(5/3-abs(z(j,k))+4*z(j,k)^2)/(5/3-abs(z(i,k))+4*z(i,k)^2);
        alpha2=(1.875-2*abs(z(j,k)))/(1.875-2*abs(z(i,k)));
        beta = (1.875-abs(z(j,k))-abs(z(Tset,k))-3*abs(z(j,k)-z(Tset,k))+8*(z(j,k)-z(Tset,k)).^2)...
            ./(1.875-abs(z(i,k))-abs(z(Tset,k))-3*abs(z(i,k)-z(Tset,k))+8*(z(i,k)-z(Tset,k)).^2);
        
        % 计算 delta
        delta = ( sum(sigma(Iset,Tset)) - sum(sigma(Jset,Tset))./beta')*(beta-1)*2;
        delta = delta+(alpha1-1)*(sum(sigma1(Jset))/alpha1-sum(sigma1(Iset)))...
            +(alpha2-1)*( (sum(sigma2(Iset))-sum(sigma2(Jset))/alpha2) );
        cc1 = 0; cc2 = 0;
        for lp = 1:KK(k)-1
            cc1 = cc1 + sum(diag(sigma(Iset,Iset),lp));
            cc2 = cc2 + sum(diag(sigma(Jset,Jset),lp));
        end
        delta = delta+2*(alpha2-1)*(cc1-cc2/alpha2);
        
        %判断是否更新
        if delta <  Disc*T-epsilon
            % 更新 D,Disc
            temp = D(i,k); D(Iset,k) = D(j,k); D(Jset,k) = temp;
            Disc = Disc + delta;
            % 更新 z,sigma1,sigma2,sigma
            temp=z(i,k);z(Iset,k)=z(j,k);z(Jset,k)=temp;            
            sigma1(Iset)=sigma1(Iset)*alpha1;
            sigma1(Jset)=sigma1(Jset)/alpha1;
            sigma2(Iset)=sigma2(Iset)*alpha2;
            sigma2(Jset)=sigma2(Jset)/alpha2;
            sigma(Tset,Iset) = sigma(Tset,Iset).*(beta*ones(1,KK(k)));
            sigma(Iset,Tset) = sigma(Tset,Iset)';
            sigma(Tset,Jset) = sigma(Tset,Jset)./(beta*ones(1,KK(k)));
            sigma(Jset,Tset) = sigma(Tset,Jset)';
            sigma(Iset,Iset) = sigma(Iset,Iset)*alpha2;
            sigma(Jset,Jset) = sigma(Jset,Jset)/alpha2;
            % 当 delta<0 说明开始下降.
            if delta < -epsilon
                % 若遇到比本次更优解，则记录之
                if Disc < MDvec(rep) - epsilon
                    MDvec(rep) = Disc;
                    % 若遇到全局更优解，则记录之
                    if MDvec(rep) < MD0-epsilon
                        MD0 = MDvec(rep);
                        if nargout > 2
                            D0 = D;
                        end
                    end
                end
            end
        end
        if writedata
            fprintf(fid,'%.8f %.8f %.8f %.8f %.8f\n',delta/Disc,T,Disc,MD0,Oiter);
        end
    end
end

if writedata
    fclose(fid);
end        

end

%{
oa = importdata('OA from web/MA.36.3.7.6.3.finney.txt');
oa = oa(:,5:10);
[N,n] = size(oa);
q = zeros(n,1);
for i = 1:n
    q(i) = length(unique(oa(:,i)));
end
[aveCD,A,BB,VEC] = aveCD_LevelPerm(oa,q);

InIter = 200;
OutIter = 200*InIter;

[MD0, MDvec, D0] = TA_uniform_FF_design_MD_abandon(oa,q,OutIter,InIter,1e-2,1e-5,10,0);
D0 = sortrows(D0);
fprintf('%.8f\n',MDvec);
%}