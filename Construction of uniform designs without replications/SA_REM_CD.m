function [Disc0, Disc_vec, ID0, D0] = SA_REM_CD(n,q,OutIter,InIter,T0,T1,beta_in,Reps,writedata)
% 20121203 
% 基于模拟退火(SA)的消除重复算法(Repetition Elimination Method)
% 求解均匀设计，迭代采用 columnwise-pairwise 策略
% INPUT:
%       n: 试验次数
%       q: s-by-1 vector, 设计各因子的水平数
%       OutIter: 外部迭代次数
%       InIter: 内部迭代次数
%       T0: 初始的温度
%       T1: 最终的温度
%       beta_in: 内部降温参数
%       Reps: 重复次数
%       writedata: 1,写数据; 0,不写(default).
% OUTPUT:
%       disc: 得到最优矩阵的 Squared Discrepancy
%       Disc_vec: Reps-by-1 vector, 记录每次重复的最优解
%       ID0: the rank vector of design D0
%       D0: 得到的最优矩阵



s = size(q,2);
if s~=1
    error('SA_REM: Input variable q must be a s-by-1 vector!\n');
end
if nargin < 8
    Reps = 1;
    writedata = 0;
elseif nargin < 9
    writedata = 0;
end
if Reps > 1
    writedata = 0;
end
% 外部降温参数计算
beta_out = (T1/T0)^(1/(OutIter-1));


s = length(q);
dd = ones(s,1);
for i=1:s-1
    dd(i)=prod(q(i+1:end));
end
epsilon = 1e-11;

N = dd(1)*q(1);

if writedata
    fid = fopen('T_Disc_Disc0.txt','w');
end

Disc0 = inf; ID0 = [];
Disc_vec = zeros(Reps,1);
for rep = 1:Reps
    % 首先生成一个随机无重复的矩阵
    [D,ID] = rand_U_Type_orth_design(q,n);
    Y = zeros(N,1); Y(ID+1) = 1;
    %生成初始的 z,sigma1,sigma
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
    % 初始化最优解信息
    % 判断当次重复的初始解是否更优
    Disc_vec(rep) = Disc;
    if Disc_vec(rep) < Disc0-epsilon
        Disc0 = Disc_vec(rep);
        if nargout > 2
            ID0=ID;
        end
    end

    Tini = T0;
    %进入模拟退火算法
    for Oiter = 1:OutIter
        T=Tini;
        ct=0;
        while ct < InIter
            ct = ct + 1;     
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
            if Y(IDi+1) || Y(IDj+1)
            %if ~isempty(find(ID==IDi,1)) || ~isempty(find(ID==IDj,1))        
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
                %if ~isempty(find(ID==IDi,1)) || ~isempty(find(ID==IDj,1))
                if Y(IDi+1) || Y(IDj+1)
                    isrepeat = 1;
                end
            end
            %产生用于计算 delta 的 alpha1,alpha2,beta, 并计算 delta
            alpha1=(1+abs(z(j,k))-2*z(j,k)^2)/(1+abs(z(i,k))-2*z(i,k)^2);
            alpha2=(1+2*abs(z(j,k)))/(1+2*abs(z(i,k)));
            beta = (1+abs(z(j,k))+abs(z(:,k))-abs(z(j,k)-z(:,k)))...
                ./(1+abs(z(i,k))+abs(z(:,k))-abs(z(i,k)-z(:,k)));
            delta = 0;
            for t = 1:n
                if t ~= i && t ~= j
                    delta = delta ...
                        + (beta(t)-1)*(sigma(i,t)-sigma(j,t)/beta(t));
                end
            end
            delta = 2*delta+(alpha1-1)*(sigma1(j)/alpha1-sigma1(i))...
                +(alpha2-1)*(sigma(i,i)-sigma(j,j)/alpha2);
        
            %判断是否更新
            if delta < -epsilon || rand() < exp(-delta/Disc/T)
                % 更新 D,ID,Disc
                Y(ID(i)+1) = 0; Y(ID(j)+1) = 0;
                Y(IDi+1) = 1; Y(IDj+1) = 1;
                temp = D(i,k); D(i,k) = D(j,k); D(j,k) = temp;
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
                % 当 delta<0 说明开始下降，为实现充分下降，设置 ct=0
                if delta < -epsilon
                    ct = 0;
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
                fprintf(fid,'%.8f %.2f %.8f %.8f %d\n',T,min(exp(-delta/Disc/T),1),Disc,Disc0,Oiter);
            end
            T = beta_in*T;
        end
        Tini = Tini*beta_out;
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
