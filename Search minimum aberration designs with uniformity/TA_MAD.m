function [aveCD0, aveCD_vec, D0] = TA_MAD(N,q,OutIter,InIter,T0,T1,Reps,filename)
%20181211 by Bochuan Jiang
%参考 Jiang and Ai (2019) Construction of uniform asymmetric designs
%采用TA算法搜索
%      q必须分组排列，即相等元素相邻排列
%INPUT:
%   N: 设计行数
%   q: n-by-1 vector, 设计各因子水平数
%   OutIter: 外部迭代次数
%   InIter: 内部迭代次数
%   T0: 初始阈值，位于(0,1)之间，建议1e-2
%   T1: 最终阈值，位于(0,1)之间，建议1e-5
%   Reps: 重复搜索次数 (default,1)
%   filename: 提供文件名的话，则将结果写入文件
%OUTPUT:
%   aveCD0: minimum average discrepancy
%   aveCD_vec: Reps重复搜索得到的 minimum average discrepancy 序列
%   D0: MAD design

N = floor(N);
if N < 1
    error('TA_MAD: Wrong N!\n');
end
[isgrouped,uq,m,group_id] = isgrouped_ByLevNum(q);
if ~isgrouped
    error('TA_MAD: q is not grouped by number of levels!\n');
end
if OutIter < 1 || InIter < 1
    error('TA_MAD: Wrong OutIter, InIter!\n');
end
if ~(0 <= T1 && T1 <= T0 && T0 < 1)
    error('TA_MAD: Wrong T0, T1!\n');
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

aveCD0 = inf;
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
for rep = 1:Reps
    % 首先生成一个随机无重复的矩阵，生成初始的 sigma
    D = rand_U_Type_design(q,N);
    sigma = zeros(N,N);
    for i = 1:N-1
        for i2 = i+1:N
            vec = abs(D(i,:)-D(i2,:))<epsilon;
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
    if aveCD_vec(rep) < aveCD0-epsilon
        aveCD0 = aveCD_vec(rep);
        if nargout > 2
            D0=D;
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
        j = randi(n);
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
                    if aveCD_vec(rep) < aveCD0-epsilon
                        aveCD0 = aveCD_vec(rep);
                        if nargout > 2
                            D0 = D;
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
%N = 16; q = [4;4;4;4;4];
%N = 36; q = [3;3;3;3;6];
N = 6; q = [2;3;3];
%N = 32; q = [4;4;4;4;4;4];
%N = 50; q = 5*ones(6,1);
InIter = 1000;
OutIter = 100*InIter;
T0 = 0.01; T1 = 1e-6;
Reps = 10;
[aveCD0, aveCD_vec, D0] = TA_MAD(N,q,OutIter,InIter,T0,T1,30,'out.txt');
aveCD = aveCD_LevelPerm(D0,q);
fprintf('%.6f  %.6f  %.6f  %d\n', aveCD0, aveCD, mean(aveCD_vec), isOA(D0,2));
aveCD_seq = importdata('out.txt');
plot(1:length(aveCD_seq),aveCD_seq);
%}



