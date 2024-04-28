function [aveDisc2, aveDisc2_vec, PI2] = TA_MAD_Combin_abandon(D0,D1,OutIter,InIter,T0,T1,flag,Reps,filename)
%20190615 by Bochuan Jiang
%参考 Jiang and Ai (2019) Construction of uniform minimum aberration designs
%采用TA局部搜索算法
%      D0固定，对 D2 进行行置换。使得 D2= [D0,D2] 的average discrepancy 最小。
%INPUT:
%   D0: N-by-n0 matrix 
%   D1: N-by-n1 matrix
%   OutIter: 外部迭代次数
%   InIter: 内部迭代次数
%   T0: 初始阈值，位于(0,1)之间，建议1e-2
%   T1: 最终阈值，位于(0,1)之间，建议1e-5
%   flag: 'CD'(default), 'WD' or 'MD'
%   Reps: 重复搜索次数 (default,1)
%   filename: 提供文件名的话，则将结果写入文件
%OUTPUT:
%   aveDisc2: minimum average discrepancy
%   aveDisc2_vec: Reps重复搜索得到的 minimum average discrepancy 序列
%   PI2: 最优的置换，D2 = [D0,D1(PI2,:)]
epsilon = 1e-10;
[N,n0] = size(D0);
[N1,n1] = size(D1);
if N ~= N1
    error('D0 and D1 have distinct numbers of rows!\n');
end
q0 = zeros(n0,1);
for i = 1:n0
    q0(i) = length(unique(D0(:,i)));
end
q1 = zeros(n1,1);
for i = 1:n1
    q1(i) = length(unique(D1(:,i)));
end
if OutIter < 1 || InIter < 1
    error('TA_MAD_Combin: Wrong OutIter, InIter!\n');
end
if ~(0 <= T1 && T1 <= T0 && T0 < 1)
    error('TA_MAD_Combin: Wrong T0, T1!\n');
end
if nargin < 7
    flag = 'CD';
end
if ~strcmp(flag,'CD') && ~strcmp(flag,'WD') && ~strcmp(flag,'MD')
    error('Wrong flag!\n');
end
if nargin < 8 || Reps < 1
    Reps = 1;
end
if nargin < 9
    iswrite = 0;
else
    iswrite = 1; % 设置输出
    fid = fopen(filename,'w');
    %记录 Reps 重复中获得最小aveCD 的搜索序列
    aveDisc_seq0 = zeros(OutIter,1);
    aveDisc_seq = zeros(OutIter,1);
    rep0 = 1; %记录最优解对应的rep
end

% 首先计算阈值变化的比率,即每 InIter 次迭代，
% 更新一次阈值 T = T*ratio
if T0 ~= 0
    ratio = (T1/T0)^(1/(ceil(OutIter/InIter)-1));
else
    ratio = 0;
end
%进入算法
q = [q0;q1];
%全局最优解 aveDisc2,PI2
aveDisc2 = aveDisc_LevelPerm([D0,D1],q,flag);
aveDisc2_vec = zeros(Reps,1);
PI2 = (1:N)';
for rep = 1:Reps    
    %产生第rep次迭代的初始解(starting design)，
    %记为循环中需要改变的量D1_loop, PI_loop, aveDisc_loop
    D1_loop = D1;
    PI_loop = (1:N)';
    aveDisc_loop = aveDisc_LevelPerm([D0,D1_loop],q,flag); 
    %第rep迭代中，局部最优的记录变量aveDisc_vec(rep),PI_opt
    aveDisc2_vec(rep) = aveDisc_loop;
    PI_opt = PI_loop;
    %如果局部最优解优于全局最优解，更新全局最优解
    if aveDisc2_vec(rep) < aveDisc2
        aveDisc2 = aveDisc2_vec(rep);
        PI2 = PI_opt;
    end
    
    %进入TA算法
    T = T0;
    for Oiter = 1:OutIter
        if mod(Oiter,InIter)==0
            T = T*ratio;
        end
        %产生D1_loop中要交换的 i,j 行
        i = randi(N);
        j = randi(N);
        while j == i
            j=randi(N);
        end
        temp = D1_loop(i,:);
        D1_loop(i,:) = D1_loop(j,:);
        D1_loop(j,:) = temp;
        
        %计算交换两行后新的 aveDisc_tilde
        aveDisc_tilde = aveDisc_LevelPerm([D0,D1_loop],q,flag);
        delta = aveDisc_tilde-aveDisc_loop;
        
        %判断是否更新
        if delta >= aveDisc_loop*T-epsilon %不更新，退回
            temp = D1_loop(i,:);
            D1_loop(i,:) = D1_loop(j,:);
            D1_loop(j,:) = temp;
        else % 更新 PI_loop,aveDisc_loop
            temp = PI_loop(i); PI_loop(i) = PI_loop(j); PI_loop(j) = temp;
            aveDisc_loop = aveDisc_tilde;
            
            % 当 delta<0 说明开始下降.
            if delta < -epsilon
                % 若遇到比本次更优解，则记录之
                if aveDisc_loop < aveDisc2_vec(rep) - epsilon
                    aveDisc2_vec(rep) = aveDisc_loop;
                    PI_opt = PI_loop;
                    % 若遇到全局更优解，则记录之
                    if aveDisc2_vec(rep) < aveDisc2-epsilon
                        aveDisc2 = aveDisc2_vec(rep);
                        PI2 = PI_opt;
                        fprintf('%d  %d  %.6f\n',rep,Oiter,aveDisc2);                       
                        if iswrite
                            rep0 = rep;
                        end
                    end
                end
            end
        end
        if iswrite 
            aveDisc_seq(Oiter) = aveDisc_loop;
        end
    end
    if iswrite && rep0 == rep
        aveDisc_seq0 = aveDisc_seq;
    end
end

if iswrite
    fprintf(fid,'%.8f\n',aveDisc_seq0);
    fclose(fid);
end

end



%{
%测试代码
InIter = 1000;
OutIter = 100*InIter;
T0 = 0.01; T1 = 1e-6;
flag = 'WD';
Reps = 1;
filename = 'out.txt';

D0 = importdata('OA from web/oa.32.9.4.2.a.txt');
D0 = importdata('OA from web/oa.18.7.3.2.txt');
D1 = D0;

[aveDisc2, aveDisc2_vec, PI2] = TA_MAD_Combin(D0,D1,OutIter,InIter,T0,T1,flag,Reps,'out.txt');
D2 = [D0,D1(PI2,:)];
aveDisc = aveDisc_LevelPerm(D2,'',flag);
fprintf('%.6f  %.6f  %.6f\n', aveDisc2, aveDisc, mean(aveDisc2_vec));
aveDisc_seq = importdata('out.txt');
plot(1:length(aveDisc_seq),aveDisc_seq);
%}



