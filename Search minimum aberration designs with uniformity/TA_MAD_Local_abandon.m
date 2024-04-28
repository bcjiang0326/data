function [aveCD1, aveCD_vec, D1] = TA_MAD_Local_abandon(D0,Jset,OutIter,InIter,T0,T1,Reps,filename)
%20181225 by Bochuan Jiang
%20190615 弃用
%参考 Jiang and Ai (2019) Construction of uniform asymmetric designs
%采用TA局部搜索算法
%      初始设计 D0 因子需按水平分组排列，即相同水平数因子相邻排列。
%      在迭代过程中仅可以改变 Jset 中元素标记的列，因此称为局部搜索。
%注：从应用角度来看，
%   (1) TA_MAD_Local 函数可以完全替代 TA_MAD_Recur 函数。比如可以在TA_MAD_Local 
%       外部构造好D0，其中，保留列是 TA_MAD_Recur 的待增加列矩阵，增加列被指定为可变列Jset。
%       这样做的好处是可以考察增加任意列情况，优于一次仅能增加一列的 TA_MAD_Recur。
%   (2) 当Jset=[1,...,n]'时，TA_MAD_Local 等同于 TA_MAD。
%INPUT:
%   D0: N-by-n matrix 一个初始设计
%   Jset: m-by-1 vector, m <= n, 元素升序排列，标记D0中可改变的列
%   OutIter: 外部迭代次数
%   InIter: 内部迭代次数
%   T0: 初始阈值，位于(0,1)之间，建议1e-2
%   T1: 最终阈值，位于(0,1)之间，建议1e-5
%   Reps: 重复搜索次数 (default,1)
%   filename: 提供文件名的话，则将结果写入文件
%OUTPUT:
%   aveCD1: minimum average discrepancy
%   aveCD_vec: Reps重复搜索得到的 minimum average discrepancy 序列
%   D1: MAD design

[N,n] = size(D0);
q = zeros(n,1);
for i = 1:n
    q(i) = length(unique(D0(:,i)));
end
[isgrouped,uq,m,group_id] = isgrouped_ByLevNum(q);
if ~isgrouped
    error('TA_MAD_Local: q is not grouped by number of levels!\n');
end
if ~issorted(Jset) || any(Jset~=floor(Jset)) || max(Jset) > n || min(Jset) < 1 ...
        || length(unique(Jset))~=length(Jset) || size(Jset,2) ~= 1
    error('TA_MAD_Local: Wrong Jset!\n');
end
if OutIter < 1 || InIter < 1
    error('TA_MAD_Local: Wrong OutIter, InIter!\n');
end
if ~(0 <= T1 && T1 <= T0 && T0 < 1)
    error('TA_MAD_Local: Wrong T0, T1!\n');
end
if nargin < 7 || Reps < 1
    Reps = 1;
end
if nargin < 8
    iswrite = 0;
else
    iswrite = 1; % 设置输出
    fid = fopen(filename,'w');
    %记录 Reps 重复中获得最小aveCD 的搜索序列
    aveCD_seq0 = zeros(OutIter,1);
    aveCD_seq = zeros(OutIter,1);
    rep0 = 1; %记录最优解对应的rep
end

n = length(q);
epsilon = 1e-11;

aveCD1 = inf;
aveCD_vec = zeros(Reps,1);

alpha1 = 2; %常数项
for k = 1:length(uq)
    if mod(uq(k),2)==0
        alpha1 = alpha1*((26*uq(k)^2+1)/24/uq(k)^2)^m(k);
    else
        alpha1 = alpha1*((13*uq(k)^2-1)/12/uq(k)^2)^m(k);
    end
end
alpha1 = (13/12)^n-alpha1;

c1 = zeros(size(uq));
for k = 1:length(uq)
    if mod(uq(k),2)==0
        c1(k) = (13*uq(k)-2)*(uq(k)-1)/12;
    else
        c1(k) = (13*uq(k)^2-2*uq(k)-3)*(uq(k)-1)/12/uq(k);
    end
end
alpha2 = prod((c1./uq./(uq-1)).^m)/N^2; %二次求和项系数

c2 = zeros(size(uq));
for k = 1:length(uq)
    if mod(uq(k),2)==0
        c2(k) = 15*uq(k)/(13*uq(k)-2);
    else
        c2(k) = (15*uq(k)^2-3)/(13*uq(k)^2-2*uq(k)-3);
    end
end



% 首先计算阈值变化的比率,即每 InIter 次迭代，
% 更新一次阈值 T = T*ratio
if T0 ~= 0
    ratio = (T1/T0)^(1/(ceil(OutIter/InIter)-1));
else
    ratio = 0;
end
%进入算法
D = D0;
for rep = 1:Reps
    % 首先对Jset标记的列进行随机化，作为第rep重复的初始矩阵
    for j = Jset'
        D(:,j) = D0(randperm(N,N),j);
    end
    % 生成初始的sigma
    sigma = zeros(N,N);
    for i = 1:N-1
        for i2 = i+1:N
            vec = abs(D(i,:)-D(i2,:)) < epsilon;
            v = zeros(length(m),1);
            %v 记录 i 行和 i2 行各个group的 Hamming distance,即相同位置的个数
            bg = 1; ed = m(1);
            v(1) = sum(vec(bg:ed));
            for k = 2:length(m)
                bg = bg+m(k-1);
                ed = bg+m(k)-1;
                v(k) = sum(vec(bg:ed));
            end
            sigma(i,i2) = alpha2*prod(c2.^v);
            sigma(i2,i) = sigma(i,i2);
        end
    end
    for i = 1:N
        sigma(i,i) = alpha2*prod(c2.^m);
    end
    
    % 初始化aveCD
    aveCD=0;
    for i = 1:N-1
        for i2 = i+1:N
            aveCD = aveCD + sigma(i,i2);
        end
    end
    aveCD = aveCD*2;
    for i = 1:N
        aveCD = aveCD + sigma(i,i);
    end
    aveCD = alpha1+aveCD;
    
    % 判断当次重复的初始解是否更优
    aveCD_vec(rep) = aveCD;
    if aveCD_vec(rep) < aveCD1-epsilon
        aveCD1 = aveCD_vec(rep);
        if nargout > 2
            D1=D;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %进入TA算法
    T = T0;
    for Oiter = 1:OutIter
        if mod(Oiter,InIter)==0
            T = T*ratio;
        end
        %产生要交换的列j及列中两个元素(i,j),(i2,j)
        j = Jset(randi(length(Jset)));
        i = randi(N);
        i2 = randi(N);
        while i2 == i || D(i,j)==D(i2,j)
            i2=randi(N);
        end
        %计算该列所在水平组
        k = group_id(j);
        
        %计算delta
        delta = 0;
        for t = 1:N
            if t ~= i && D(t,j)==D(i,j)
                delta = delta + sigma(i2,t)-sigma(i,t)/c2(k);
            elseif t~=i2 && D(t,j)==D(i2,j)
                delta = delta + sigma(i,t)-sigma(i2,t)/c2(k);
            end
        end
        delta = delta*2*(c2(k)-1);
        
        %判断是否更新
        if delta <  aveCD*T-epsilon
            % 更新 sigma
            for t = 1:N
                if t ~= i && D(t,j) == D(i,j)
                    sigma(i,t) = sigma(i,t)/c2(k);
                    sigma(i2,t) = sigma(i2,t)*c2(k);
                    sigma(t,i)=sigma(i,t);
                    sigma(t,i2)=sigma(i2,t);
                elseif t ~= i2 && D(t,j) == D(i2,j)
                    sigma(i,t) = sigma(i,t)*c2(k);
                    sigma(i2,t) = sigma(i2,t)/c2(k);
                    sigma(t,i)=sigma(i,t);
                    sigma(t,i2)=sigma(i2,t);
                end
            end
            % 更新 D,aveCD
            temp=D(i,j);D(i,j)=D(i2,j);D(i2,j)=temp;
            aveCD = aveCD + delta;
            % 当 delta<0 说明开始下降.
            if delta < -epsilon
                % 若遇到比本次更优解，则记录之
                if aveCD < aveCD_vec(rep) - epsilon
                    aveCD_vec(rep) = aveCD;
                    % 若遇到全局更优解，则记录之
                    if aveCD_vec(rep) < aveCD1-epsilon
                        aveCD1 = aveCD_vec(rep);
                        if nargout > 2
                            D1 = D;
                        end
                        if iswrite
                            rep0 = rep;
                        end
                    end
                end
            end
        end
        if iswrite 
            aveCD_seq(Oiter) = aveCD;
        end
    end
    if iswrite && rep0 == rep
        aveCD_seq0 = aveCD_seq;
    end
end

if iswrite
    fprintf(fid,'%.8f\n',aveCD_seq0);
    fclose(fid);
end

end



%{
%测试代码

k = 2; m = 2; p = 2; 
s = p^m; 
prim_poly = gfprimfd(m,'min',p);
[D0,GenMat,OA_polys,GM_polys] = RH_OA_pow(k,m,p,prim_poly);
[N,n] = size(D0);
q = s*ones(n,1);
Jset = [n-1;n];

InIter = 1000;
OutIter = 100*InIter;
T0 = 0.01; T1 = 1e-6;
Reps = 10;
filename = 'out.txt';

[aveCD1, aveCD_vec, D1] = TA_MAD_Local(D0,Jset,OutIter,InIter,T0,T1,Reps,'out.txt');
aveCD = aveCD_LevelPerm(D1,q);
fprintf('%.6f  %.6f  %.6f  %d\n', aveCD1, aveCD, mean(aveCD_vec), isOA(D1,2));
aveCD_seq = importdata('out.txt');
plot(1:length(aveCD_seq),aveCD_seq);
%}



