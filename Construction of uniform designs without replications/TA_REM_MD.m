function [Disc0, Disc_vec, ID0, D0] = TA_REM_MD(n,q,OutIter,InIter,T0,T1,Reps,writedata)
% 20121203 
% 基于阈值接受(TA)的消除重复算法(Repetition Elimination Method)
% The threshold is lowered gemetrically after every InIter itertations
% 求解均匀设计，迭代采用 columnwise-pairwise 策略
% INPUT:
%       n: 试验次数
%       q: s-by-1 vector, 设计各因子的水平数
%       OutIter: 外部迭代次数
%       InIter: 内部迭代次数
%       T0: 初始阈值，位于(0,1)之间，建议1e-2
%       T1: 最终阈值，位于(0,1)之间，建议1e-5
%       Reps: 重复次数 (default,1)
%       writedata: 是否将过程数据写入文件 -- 1,写入; 0(default)
% OUTPUT:
%       Disc0: 得到最优矩阵的 Squared Discrepancy
%       Disc_vec: Reps-by-1 vector, 记录每次重复的最优解
%       ID0: the rank vector of design D0
%       D0: 得到的最优矩阵
s = size(q,2);
if s~=1
    error('TA_REM_MD: Input variable q must be a s-by-1 vector!\n');
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
% 首先计算阈值变化的比率,即每 InIter 次迭代，
% 更新一次阈值 T = T*ratio     
ratio = (T1/T0)^(1/(ceil(OutIter/InIter)-1));

s = length(q);
dd = ones(s,1);
for i=1:s-1
    dd(i)=prod(q(i+1:end));
end
epsilon = 1e-11;

% 设置输出
if writedata
    fid = fopen('T_Disc_Disc0.txt','w');
end

Disc0 = inf; ID0 = [];
Disc_vec = zeros(Reps,1);

for rep = 1:Reps
    % 首先生成一个随机无重复的矩阵，生成初始的 z,sigma1,sigma
    [D,ID] = rand_U_Type_orth_design(q,n);
    z = zeros(n,s);
    sigma1 = zeros(n,1);
    sigma = zeros(n,n);
    for i = 1:n
        z(i,:) = 0.5*( (D(i,:)+0.5)./q'-0.5 );
    end
    for i = 1:n-1
        for j = i+1:n
            sigma(i,j) = 1/n^2;
            for k = 1:s
                sigma(i,j)=sigma(i,j)*(1+abs(z(i,k))+abs(z(j,k))-abs(z(i,k)-z(j,k)));
            end
            sigma(j,i)=sigma(i,j);
        end 
    end
    for i = 1:n
        sigma1(i) = 2/n;
        for k = 1:s
            sigma1(i) = sigma1(i)*(1+abs(z(i,k))-2*z(i,k)^2);
        end
        sigma(i,i) = prod(1+2*abs(z(i,:)))/n^2;
    end
    
    % 初始化Disc
    Disc=0;
    for i = 1:n-1
        for j = i+1:n
            Disc = Disc + sigma(i,j);
        end
    end
    Disc=(13/12)^s+2*Disc+sum(diag(sigma))-sum(sigma1);

    % 判断当次重复的初始解是否更优
    Disc_vec(rep) = Disc;
    if Disc_vec(rep) < Disc0-epsilon
        Disc0 = Disc_vec(rep);
        if nargout > 2
            ID0=ID;
        end
    end

    %进入算法
    T = T0;
    for Oiter = 1:OutIter
        if mod(Oiter,InIter)==0
            T = T*ratio;
        end
        %产生要交换的列及列中两个元素    
        isrepeat = 0;
        k=randi(s);
        i=randi(n);
        j=randi(n);
        while j==i || D(i,k)==D(j,k)    
            j=randi(n);
        end
        IDi = ID(i) + (D(j,k)-D(i,k))*dd(k); %此处记录下此二数，
        IDj = ID(j) + (D(i,k)-D(j,k))*dd(k); %方便后面更新 ID
        if ~isempty(find(ID==IDi,1)) || ~isempty(find(ID==IDj,1))        
            isrepeat = 1;        
        end        
        %此时首先考虑交换是否会带来重复行
        %若出现重复行，则重新进行交换
        while isrepeat
            isrepeat = 0;
            k=randi(s);
            i=randi(n);
            j=randi(n);
            while j==i || D(i,k)==D(j,k)
                j=randi(n);
            end
            %检验是否可能出现重复
            IDi = ID(i) + (D(j,k)-D(i,k))*dd(k);
            IDj = ID(j) + (D(i,k)-D(j,k))*dd(k);
            if ~isempty(find(ID==IDi,1)) || ~isempty(find(ID==IDj,1))
                isrepeat = 1;
            end
        end
    
        %产生用于计算 delta 的 alpha1,alpha2,beta
        alpha1=(1+abs(z(j,k))-2*z(j,k)^2)/(1+abs(z(i,k))-2*z(i,k)^2);
        alpha2=(1+2*abs(z(j,k)))/(1+2*abs(z(i,k)));
        beta = (1+abs(z(j,k))+abs(z(:,k))-abs(z(j,k)-z(:,k)))...
            ./(1+abs(z(i,k))+abs(z(:,k))-abs(z(i,k)-z(:,k)));
    
        %计算delta          
        delta = 0;
        for t = 1:n
            if t ~= i && t ~= j
                delta = delta ...
                    +(beta(t)-1)*(sigma(i,t)-sigma(j,t)/beta(t));
            end
        end
        delta = 2*delta+(alpha1-1)*(sigma1(j)/alpha1-sigma1(i))...
            +(alpha2-1)*(sigma(i,i)-sigma(j,j)/alpha2);
        
        %判断是否更新
        if delta <  Disc*T-epsilon
            % 更新 D,ID,Disc
            temp=D(i,k);D(i,k)=D(j,k);D(j,k)=temp;
            ID(i) = IDi; ID(j) = IDj;
            Disc = Disc + delta;
            % 更新 z,sigma1,diag(sigma),sigma
            temp=z(i,k);z(i,k)=z(j,k);z(j,k)=temp;            
            sigma1(i)=sigma1(i)*alpha1;
            sigma1(j)=sigma1(j)/alpha1;
            sigma(i,i)=sigma(i,i)*alpha2;
            sigma(j,j)=sigma(j,j)/alpha2;
            for t=1:n
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
                            ID0 = ID;
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
if nargout > 2
    ID0 = sort(ID0);
end
if nargout > 3
    D0 = Id2Design(q,ID0);
end
end


