function [WD0, WDvec, D0] = TA_uniform_FF_design_inverseWD(D,AllIter,InIter,T0,T1,Reps,writedata)
% 20200213 求WD最大的设计，其余参考函数 TA_uniform_FF_design(...)      .
% INPUT:
%       D: N-by-n balanced design
%       AllIter: 总迭代次数=OutIter*InIter
%       InIter: 内部迭代次数
%       T0: 初始阈值，位于(0,1)之间，建议1e-2
%       T1: 最终阈值，位于(0,1)之间，建议1e-6
%       Reps: 重复次数 (default,1)
%       writedata: 是否将过程数据写入文件。 1,写入; 0(default) 不写
% OUTPUT:
%       WD0: 得到最优矩阵的 Squared WD
%       WDvec: Reps-by-1 vector, 记录每次重复的最优解
%       D0: 得到的最优矩阵

[N,n] = size(D);
[isbal,q] = isBalanced(D);
if ~isbal
    error('D is not balanced!\n');
end
if nargin < 6
    Reps = 1;
    writedata = 0;
elseif nargin < 7
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
epsilon = 1e-12;

% 首先计算阈值变化的比率,即每 InIter 次迭代，
% 更新一次阈值 T = T*ratio
ratio = (T1/T0)^(1/(ceil(AllIter/InIter)-1));

KK = N./q; % 记录各列中，每个水平重复次数

D0 = D;
WD0 = WD2_value(D0,q);
WDvec = zeros(Reps,1);

for rep = 1:Reps
    % level_id{j}记录第j列中各元素出现的位置
    level_id = cell(1,n);
    for j = 1:n
        [~,id] = sort(D(:,j));
        level_id{j} = reshape(id,[N/q(j),q(j)]);
    end
    % 首先生成初始的 sigma
    sigma = zeros(N,N);
    for i = 1:N-1
        for j = i+1:N
            avec = abs(D(i,:)-D(j,:))./q';
            sigma(i,j) = 1/N^2*prod(1.5-avec.*(1-avec));
            sigma(j,i)=sigma(i,j);
        end
    end
    
    % 初始化Disc
    Disc=0;
    for i = 1:N-1
        for j = i+1:N
            Disc = Disc + sigma(i,j);
        end
    end
    Disc=-(4/3)^n + 1.5^n/N + 2*Disc;
    
    % 判断当次重复的初始解是否更优
    WDvec(rep) = Disc;
    if WDvec(rep) > WD0+epsilon
        WD0 = WDvec(rep);
        if nargout > 2
            D0=D;
        end
    end
    
    %进入算法
    T = T0;
    for Oiter = 1:AllIter
        if mod(Oiter,InIter)==0
            T = T*ratio;
        end
        %产生要交换的列
        k=randi(n);
        %k = n;
        %产生要交换的两个水平lv1,lv2
        lv1=randi(q(k))-1;
        lv2=randi(q(k))-1;
        while lv2==lv1
            lv2=randi(q(k))-1;
        end
        if lv1 > lv2
            temp = lv1; lv1 = lv2; lv2 = temp;
        end
        %置换2水平因子、3水平因子、或者4水平因子的0和2、或者4水平因子的1和3，均不改变WD值.
        if q(k)==2 || q(k)==3 || ( q(k) == 4 && lv1 == 0 && lv2 == 2 ) || ( q(k) == 4 && lv1 == 1 && lv2 == 3 )
            if writedata
                fprintf(fid,'%.8f %.8f %.8f %.8f %.8f\n',0,T,Disc,WD0,Oiter);
            end
            continue;
        end
        Iset = level_id{k}(:,lv1+1);
        Jset = level_id{k}(:,lv2+1);
        Tset = reshape(level_id{k}(:,[1:lv1,lv1+2:lv2,lv2+2:end]),[N-2*KK(k),1]);
        
        %产生用于计算 delta 的 beta
        ivec = abs(lv1-D(Tset,k))/q(k);
        jvec = abs(lv2-D(Tset,k))/q(k);
        beta = (1.5 - jvec.*(1-jvec))./(1.5-ivec.*(1-ivec));
        
        % 计算 delta
        delta = ( sum(sigma(Iset,Tset)) - sum(sigma(Jset,Tset))./beta')*(beta-1)*2;
        
        %判断是否更新
        uu = rand();
        if delta >  -Disc*T*(1+epsilon) || uu < 1e-3
            % 更新 D,Disc,level_id
            D(Iset,k) = lv2; D(Jset,k) = lv1;
            Disc = Disc + delta;
            level_id{k}(:,lv1+1) = Jset;
            level_id{k}(:,lv2+1) = Iset;
            % 更新 sigma
            sigma(Tset,Iset) = sigma(Tset,Iset).*(beta*ones(1,KK(k)));
            sigma(Iset,Tset) = sigma(Tset,Iset)';
            sigma(Tset,Jset) = sigma(Tset,Jset)./(beta*ones(1,KK(k)));
            sigma(Jset,Tset) = sigma(Tset,Jset)';
            % 当 delta>0 说明开始上升.
            if delta > Disc*epsilon
                % 若遇到比本次更大解，则记录之
                if Disc > WDvec(rep)*(1 + epsilon)
                    WDvec(rep) = Disc;
                    % 若遇到全局更优解，则记录之
                    if WDvec(rep) > WD0*(1 + epsilon)
                        WD0 = WDvec(rep);
                        D0 = D;
                        fprintf('inv LP: N=%d, n=%d, rep=%d, iter=%d: %.6f\n',N,n,rep,Oiter,WD0);
                    end
                end
            end
        end
        if writedata
            fprintf(fid,'%.8f %.8f %.8f %.8f %.8f\n',delta/Disc,T,Disc,WD0,Oiter);
        end
    end
    %fprintf('LP: N = %d, n = %d, rep = %d: ******** %.6f, %.6f\n',N,n,rep,WDvec(rep),WD0);
end

if writedata
    fclose(fid);
end

end

%{
oa = importdata('OA from web/MA.36.3.7.6.3.finney.txt');
oa = oa(:,5:10);
aveWD = aveDisc_LevelPerm(oa,'','WD');

InIter = 100;
OutIter = 100*InIter;

[WD0, WDvec, D0] = TA_uniform_FF_design_inverseWD(oa,OutIter,InIter,1e-3,1e-6,1,1);
D0 = sortrows(D0);
fprintf('%.8f\n',WDvec);
fprintf('%.8f\n',WD2_value(D0));
out = importdata('T_Disc_Disc.txt');
plot(out(:,5),out(:,3));
%}