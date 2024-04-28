function [ptn0, ptn_vec, D0] = TA_proj_CD0(dim, N,q,OutIter,InIter,T0,T1,Reps,writedata)
% 20131001 最初版本，dim 只能为标量
% 基于 TA 算法求解特定维度（dim）下投影均匀性
% 求解均匀设计，迭代采用 columnwise-pairwise 策略
% 背景：Hickernell and Liu (2002) 提出 projection discrepancy pattern,
%       本程序计算特别维度投影均匀性最优的设计。程序框架采用 TA 算法。
%       .
% INPUT:
%       dim：scalar 投影维度
%       N: run-size
%       q: n-by-1 vector, 设计各因子的水平数
%       OutIter: 外部迭代次数
%       InIter: 内部迭代次数
%       T0: 初始阈值，位于(0,1)之间，建议1e-2
%       T1: 最终阈值，位于(0,1)之间，建议1e-6
%       Reps: 重复次数 (default,1)
%       writedata: 是否将过程数据写入文件。 1,写入; 0(default) 不写
% OUTPUT:
%       ptn0: 得到最优矩阵的 Squared Discrepancy
%       ptn_vec: Reps-by-1 vector, 记录每次重复的最优解
%       D0: 得到的最优矩阵

if min(size(q)) ~= 1 || N ~= floor(N) || any( N./q ~= floor(N./q) )
    error(' Input variable q must be a s-by-1 vector!\n');
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
n = length(q);

% 设置文件输出
if writedata
    fid = fopen('T_Disc_Disc0.txt','w');
end

% 设置归零参数
epsilon = 1e-12;

% 首先计算阈值变化的比率,即每 InIter 次迭代，
% 更新一次阈值 T = T*ratio     
ratio = (T1/T0)^(1/(ceil(OutIter/InIter)-1));


ptn0 = inf; 
ptn_vec = zeros(Reps,1);

for rep = 1:Reps
    % 首先生成一个随机无重复的矩阵
    D = rand_U_Type_design(q,N);
    DD = (D+0.5)./(2*ones(N,1)*q')-0.25;
    aD = abs(DD);
    
    % 初始化 PDisc
    ptn = CD2_pattern(D,q,dim);

    % 判断当次重复的初始解是否更优
    ptn_vec(rep) = ptn;
    if ptn < ptn0-epsilon
        ptn0 = ptn;
        if nargout > 2
            D0=D;
        end
    end
    
    %进入算法
    T = T0;
    for Oiter = 1:OutIter
        if mod(Oiter,InIter)==0
            T = T*ratio;
        end
        
        %产生要交换的列及列中两个元素    
        k0 = randi( n ); % 随机选第 k0 列
        i0 = randi( N );
        j0 = randi( N );
        while j0==i0 || D(i0,k0) == D(j0,k0)    
            j0=randi( N );
        end
        
        %i0 = 2; j0 = 4; k0 = 1;
        
   
        %计算delta(kk) = delta1 + delta21 + delta22
        delta1 = 0; delta21 = 0; delta22 = zeros(1,N);
        id  = 1:dim-1;
        while id(1)~=-1
            id0 = id;
            id0(id0>=k0) = id0(id0>=k0)+1;
            delta1 = delta1+prod(aD(i0,id0)-2*DD(i0,id0).^2)-prod(aD(j0,id0)-2*DD(j0,id0).^2);
            delta21 = delta21+prod(aD(i0,id0))-prod(aD(j0,id0));
            for t = 1:N
                if t~=i0 && t~=j0
                    delta22(t) = delta22(t)+prod( aD(i0,id0)+aD(t,id0)-abs(DD(i0,id0)-DD(t,id0)) ) ...
                        -prod( aD(j0,id0)+aD(t,id0)-abs(DD(j0,id0)-DD(t,id0)) );
                end
            end
            id = nchoosek_next(id,n-1,dim-1);
        end
        delta1 = delta1*(-2/N)*(aD(j0,k0)-aD(i0,k0))*(1-2*(aD(j0,k0)+aD(i0,k0)));
        delta21 = delta21*2^dim/N^2*(aD(j0,k0)-aD(i0,k0));
        delta22 = delta22*(aD(j0,k0)-aD(i0,k0)-abs(DD(j0,k0)-DD(:,k0))+abs(DD(i0,k0)-DD(:,k0)))*2/N^2;
        delta = delta1+delta21+delta22;
        
        %判断是否更新
        %if delta <  PDisc*T-epsilon
        if  (rand()<0.001 && T < 1e-4) ||  delta<  ptn*T-epsilon
            % 更新 D,DD,PDisc
            temp=D(i0,k0);D(i0,k0)=D(j0,k0);D(j0,k0)=temp;
            temp=DD(i0,k0);DD(i0,k0)=DD(j0,k0);DD(j0,k0)=temp;
            temp=aD(i0,k0); aD(i0,k0)=aD(j0,k0); aD(j0,k0)=temp;
            ptn = ptn + delta;
            % 当 delta<0 说明开始下降.
            if delta < -epsilon
                % 若遇到比本次更优解，则记录之
                if ptn < ptn_vec(rep)-epsilon
                    ptn_vec(rep) = ptn;
                    % 若遇到全局更优解，则记录之
                    if ptn_vec(rep) < ptn0
                        ptn0 = ptn_vec(rep);
                        if nargout > 2
                            D0 = D;
                        end
                    end
                end
            end
        end
        if writedata
            fprintf(fid,'%.8f %.8f %.8f %.8f %.8f\n',delta/ptn,T,ptn,ptn0,Oiter);
        end
    end
end

if writedata
    fclose(fid);
end
end

%{
% 测试函数
D=[2 4 2 4; 4 3 2 1; 3 4 4 2; 2 1 3 1; 4 2 3 4; 3 1 1 3; 1 2 1 2; 1 3 4 3];
D = sortrows(D-1); % Uniform (8,4^4)-design

s = max(max(D))+1;
[N,n] = size(D);
q = s*ones(n,1);


InIter = 200;
OutIter = 200*InIter;
Reps = 1;
dims = [2,3];
[ptn0, ptn_vec, D0] = TA_proj_CD(dims,N,q,OutIter,InIter,1e-2,1e-6,Reps,0);
fprintf('%.4e ',ptn0);fprintf('\n');
%}