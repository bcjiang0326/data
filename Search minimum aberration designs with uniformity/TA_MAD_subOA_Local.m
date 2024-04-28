function [aveCD1, aveCD_vec,Jin1,D1] = TA_MAD_subOA_Local(D0,q1,JinFix,OutIter,InIter,T0,T1,Reps,filename)
%20181225 by Bochuan Jiang
%采用TA局部搜索算法
%      D0为一个已知正交表，搜索其水平为q1序列的 MAD 子表（sub-OA）
%      初始设计 D0 因子需按水平分组排列，即相同水平数因子相邻排列。
%      JinFix为必须要放入sub-OA中的列标号集合，在迭代过程中，JinFix
%      中的列不可被移出，因此称为局部搜索。
%注意：（1）当T0=0时，TA算法变成了Greedy算法
%     （2）本函数是从正交表中搜索MAD子表的一般性函数，包含类似功能的其它函数。
%INPUT:
%   D0: N-by-n matrix 一个初始设计
%   q1: n1-by-1 vector, sub-OA 的各个水平
%   JinFix: m-by-1 vector, m <= n, 元素升序排列，必须要放入 subOA 的 D0 的列标号
%   OutIter: 外部迭代次数
%   InIter: 内部迭代次数
%   T0: 初始阈值，位于(0,1)之间，建议1e-2
%   T1: 最终阈值，位于(0,1)之间，建议1e-5
%   Reps: 重复搜索次数 (default,1)
%   filename: 提供文件名的话，则将结果写入文件
%OUTPUT:
%   aveCD1: minimum average discrepancy
%   aveCD_vec: Reps重复搜索得到的 minimum average discrepancy 序列
%   Jin1: D0(:,Jin1) = D1
%   D1: sub-OA

[N,n] = size(D0);
q0 = zeros(n,1);
for i = 1:n
    q0(i) = length(unique(D0(:,i)));
end
[isgrouped,uq,m,group_id] = isgrouped_ByLevNum(q0);
if ~isgrouped
    error('TA_MAD_subOA_Local: q is not grouped by number of levels!\n');
end
if ~issorted(JinFix) || any(JinFix~=floor(JinFix)) || max(JinFix) > n || min(JinFix) < 1 ...
        || length(unique(JinFix))~=length(JinFix) || size(JinFix,2) ~= 1
    error('TA_MAD_subOA_Local: Wrong Jin!\n');
end
m1 = zeros(size(m));
for k = 1:length(uq)
    m1(k) = sum(q1==uq(k));
end
if any(m1>m) %相当于比较q1和q
    error('TA_MAD_subOA_Local: Wrong q1!\n');
end
m_Jin = zeros(size(m));
for k = 1:length(JinFix)
    g = group_id(JinFix(k));
    m_Jin(g) = m_Jin(g)+1;
end
if any(m_Jin>m1)
    error('TA_MAD_subOA_Local: q1 and Jin do not match!\n');
end
if OutIter < 1 || InIter < 1
    error('TA_MAD_subOA_Local: Wrong OutIter, InIter!\n');
end
if ~(0 <= T1 && T1 <= T0 && T0 < 1)
    error('TA_MAD_subOA_Local: Wrong T0, T1!\n');
end
if nargin < 8 || Reps < 1
    Reps = 1;
end
if length(JinFix) == length(q1)
    %这种情况下，JinFix对应的就是最终结果
    Jin1 = JinFix;
    D1 = D0(:,JinFix);
    aveCD1 = aveCD_LevelPerm(D1);
    aveCD_vec = aveCD1*ones(Reps,1);
    return;
end
if nargin < 9
    iswrite = 0;
else
    iswrite = 1; % 设置输出
    fid = fopen(filename,'w');
    %记录 Reps 重复中获得最小aveCD 的搜索序列
    aveCD_seq0 = zeros(OutIter,1);
    aveCD_seq = zeros(OutIter,1);
    rep0 = 1; %记录最优解对应的rep
end

%设置精度和输出变量初值
epsilon = 1e-11;
aveCD1 = inf;
aveCD_vec = zeros(Reps,1);
Jin1 = []; % 用Jin1 来记录全局（所有重复下）最优sub-OA

alpha1 = 2; %常数项
for k = 1:length(uq)
    if mod(uq(k),2)==0
        alpha1 = alpha1*((26*uq(k)^2+1)/24/uq(k)^2)^m1(k);
    else
        alpha1 = alpha1*((13*uq(k)^2-1)/12/uq(k)^2)^m1(k);
    end
end
alpha1 = (13/12)^(length(q1))-alpha1;

c1 = zeros(size(uq));
for k = 1:length(uq)
    if mod(uq(k),2)==0
        c1(k) = (13*uq(k)-2)*(uq(k)-1)/12;
    else
        c1(k) = (13*uq(k)^2-2*uq(k)-3)*(uq(k)-1)/12/uq(k);
    end
end
alpha2 = prod((c1./uq./(uq-1)).^m1)/N^2; %二次求和项系数

c2 = zeros(size(uq));
for k = 1:length(uq)
    if mod(uq(k),2)==0
        c2(k) = 15*uq(k)/(13*uq(k)-2);
    else
        c2(k) = (15*uq(k)^2-3)/(13*uq(k)^2-2*uq(k)-3);
    end
end

%计算D0中任意两列的Iset:
%Iset{j1,j2}={(i1,i2):i1>i2 and D0(i1,j1)=D0(i2,j1) and D0(i1,j2)!=D0(i2,j2)}
%代表j1列中相等的行对同时j2列中不等的行对。
Iset = cell(n,n);
for j1 = 2:n
    for j2 = 1:j1-1
        if group_id(j1)~=group_id(j2)
            %非同一水平组变量，跳过
            continue;
        end
        Iset{j1,j2} = zeros(0,2);
        Iset{j2,j1} = zeros(0,2);
        %如果D0是强度大于等于2的正交表，则Iset的长度为 (N/s)^2(s-1)/2
        for i1 = 2:N
            for i2 = 1:i1-1
                if D0(i1,j1)==D0(i2,j1) && D0(i1,j2)~=D0(i2,j2)
                    Iset{j1,j2} = cat(1,Iset{j1,j2},[i1,i2]);
                elseif D0(i1,j2)==D0(i2,j2) && D0(i1,j1)~=D0(i2,j1)
                    Iset{j2,j1} = cat(1,Iset{j2,j1},[i1,i2]);
                end
            end
        end
    end
end


% 首先计算阈值变化的比率,即每 InIter 次迭代，
% 更新一次阈值 T = T*ratio
if T0 ~= 0
    ratio = (T1/T0)^(1/(ceil(OutIter/InIter)-1));
else
    ratio = 0;
end

%计算每个group起始列的上一列标号
bg = zeros(length(m),1);
bg(2:end) = cumsum(m(1:end-1));

%进入算法
for rep = 1:Reps
    % 首先对Jset标记的列进行随机化，作为第rep重复的初始矩阵
    Jin = JinFix;
    m_left = m1-m_Jin;
    for k = 1:length(m1)
        while m_left(k) > 0
            j = bg(k)+randi(m(k));
            while any(j==Jin)
                j = bg(k)+randi(m(k));
            end
            Jin = cat(1,Jin,j);
            m_left(k) = m_left(k)-1;
        end
    end
    Jin = sort(Jin);

    % 生成初始的sigma
    sigma = zeros(N,N);
    for i = 1:N-1
        for i2 = i+1:N
            vec = abs(D0(i,Jin)-D0(i2,Jin)) < epsilon;
            v = zeros(length(m1),1);
            %v 记录 i 行和 i2 行各个group的 Hamming distance,即相同位置的个数
            bgg = 1; edd = m1(1);
            v(1) = sum(vec(bgg:edd));
            for k = 2:length(m1)
                bgg = bgg+m1(k-1);
                edd = bgg+m1(k)-1;
                v(k) = sum(vec(bgg:edd));
            end
            sigma(i,i2) = alpha2*prod(c2.^v);
            sigma(i2,i) = sigma(i,i2);
        end
    end
    for i = 1:N
        sigma(i,i) = alpha2*prod(c2.^m1);
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
        Jin1 = Jin;
    end
    
    %进入TA算法
    T = T0;
    for Oiter = 1:OutIter
        if mod(Oiter,InIter)==0
            T = T*ratio;
        end
        %产生要移出的列j1
        j1_id = randi(length(Jin));
        j1 = Jin(j1_id);
        while any(j1==JinFix)
            j1_id = randi(length(Jin));
            j1 = Jin(j1_id);
        end
        k = group_id(j1); %该列所在水平组
        %产生要移入的列j2
        j2 = bg(k)+randi(m(k));
        while any(j2==Jin)
            j2 = bg(k)+randi(m(k));
        end
        
        %计算delta
        delta1 = 0;
        for t = 1:size(Iset{j2,j1},1)
            i1 = Iset{j2,j1}(t,1);
            i2 = Iset{j2,j1}(t,2);
            delta1 = delta1+sigma(i1,i2);
        end
        delta2 = 0;
        for t = 1:size(Iset{j1,j2},1)
            i1 = Iset{j1,j2}(t,1);
            i2 = Iset{j1,j2}(t,2);
            delta2 = delta2+sigma(i1,i2);
        end
        delta = 2*(c2(k)-1)*(delta1-delta2/c2(k));
        
        %判断是否更新
        if delta <  aveCD*T-epsilon
            % 更新 sigma
            for t = 1:size(Iset{j2,j1},1)
                i1 = Iset{j2,j1}(t,1);
                i2 = Iset{j2,j1}(t,2);
                sigma(i1,i2) = sigma(i1,i2)*c2(k);
                sigma(i2,i1) = sigma(i1,i2);
            end
            for t = 1:size(Iset{j1,j2},1)
                i1 = Iset{j1,j2}(t,1);
                i2 = Iset{j1,j2}(t,2);
                sigma(i1,i2) = sigma(i1,i2)/c2(k);
                sigma(i2,i1) = sigma(i1,i2);
            end
            % 更新Jin,aveCD
            Jin(j1_id) = j2; Jin = sort(Jin);
            aveCD = aveCD + delta;
            % 当 delta<0 说明开始下降.
            if delta < -epsilon
                % 若遇到比本次更优解，则记录之
                if aveCD < aveCD_vec(rep) - epsilon
                    aveCD_vec(rep) = aveCD;
                    % 若遇到全局更优解，则记录之
                    if aveCD_vec(rep) < aveCD1-epsilon
                        aveCD1 = aveCD_vec(rep);
                        Jin1 = Jin;
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
if nargout > 3
    D1 = D0(:,Jin1);
end

end



%{
%测试代码
InIter = 100;
OutIter = 100*InIter;
T0 = 1e-2; T1 = 1e-5;
Reps = 1;
filename = 'out.txt';

%混水平情形
D0 = importdata('OA from web/MA.36.2.11.3.12.txt');
q1 = [2;3;2;3;2];
JinFix = [1;12];
[aveCD1,aveCD_vec,Jin1,D1] = TA_MAD_subOA_Local(D0,q1,JinFix,OutIter,InIter,T0,T1,Reps,'out.txt');
data = importdata('out.txt');
figure
plot(1:size(data,1),data(:,1));

%等水平情形
k = 3; s = 3; 
[D0,GenMat] = RH_OA_prim(k,s);
JinFix = [1];
%JinFix = [1;2;5;14];
%JinFix = [1;2;5;8;14;17;20;21];
n1 = 3;
q1 = s*ones(n1,1);
[aveCD1,aveCD_vec,Jin1,D1] = TA_MAD_subOA_Local(D0,q1,JinFix,OutIter,InIter,T0,T1,Reps,'out.txt');
[aveCD_,A_] = aveCD_LevelPerm(D1);
data = importdata('out.txt');
figure
plot(1:size(data,1),data(:,1));
%}



