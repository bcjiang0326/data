function [x wd2 info] = UD_UTYPE(Fun,q0,beq,x0,upperbnd,lowerbnd,opt0)
% 2012/4/24 2012/4/25
% 剪枝技术：1.Convex optimization；2.t；3.U-type design
% 顺序分裂，因此在递归过程中不需要向量 vec。取而代之的是标量 j
% Input
%       Fun: function_handle, 计算 H 和 f：0.5*x'*H*x+f'*x, 不同
%            的均匀性准则对应不同的 Fun.
%       q0: 各个因子的水平数 s-by-1 向量
%       x0: 必须是一个初始的可行解并且 bound = 0.5*x'*H*x
%       beq: 约束条件 sum(x) = beq 并且 beq < m
%       upperbnd: 到目前为止已找到的最优值。 default = inf
%       lowerbnd: 先验下界，未必能够达到。  default = -inf
%       opt0: 
%           opt0(1):是否采用凸优化剪枝。是 = 1,default 0
%           opt0(2):是否采用 wd2Ls 中的 neighborhood 技术 进行底层运算。是 = 1，default = 0
%           opt0(3):采用 t 剪枝的阙值，范围(0,1]。default = 1 表示不采用。
% Output
%       x: 搜索获得的最优解
%       info: 存储各项信息的数组
%           info(1): 在第一处更新 upperbnd 的次数
%           info(2): 在第二处更新 upperbnd 的次数
%           info(3): 递归调用次数
%           info(4): upperbnd 初始值
%           info(5): upperbnd 终止值 
%           info(6): upperbnd 相对于 lowerbnd 的误差
%                   i.e. (info(5)-lowerbnd)/lowerbnd.
%                   取-1, 若 lowerbnd <= 0.    
%           info(7): 访问的叶节点个数



narginchk(3, 7);%, nargin, 'struct'));
nargoutchk(0, 3);%, nargout, 'struct'));

% 全局变量声明 部分1
% H: 输入矩阵，在上级目录中计算
% q: 各个因子水平组合数 s-by-1 
% TolXInt: 判断为整数的阙值，在凸优化后用
% options: 凸优化中的参数
% xx: 初始点，并存储最终结果
% epsilon: 本算法的精度
% n: 总的试验次数
% m: 因子水平组合总数
% s: 因子个数
global H q TolXInt options xx epsilon info1 n m s opt;

% 全局变量声明 部分2
% TriesPerLel：记录每个因子的每个水平上已经试验的次数，sum(q)-by-1 向量
% NoTriesPerLel：记录每个因子的每个水平上安排的0的个数，sum(q)-by-1 向量
% URestrict：s-by-4 矩阵，
%       第一列：各因子平均每水平试验次数（向下取整）i.e.分配1的个数 floor(beq./q)
%       第二列：各因子试验允许超出平均的水平数上限 mod(beq,q);
%       第三列：各因子平均每水平不实验的次数（向下取整）i.e.分配0的个数 floor((m-beq)./q)
%       第四列：各因子不试验次数允许超出平均的水平数上限  mod(m-beq,q);
% Staff: 标尺，各个因子在 TriesPerLel 中分界的位置，起始位置前一点
% NGetRestr1：记录当前各个因子中已经超出平均的水平数
% NGetRestr0: 记录各个因子中不试验次数已经超出平均的水平数
% 不试验次数 —— 0 的个数；试验次数 —— 1 的个数。
global TriesPerLel NoTriesPerLel URestrict Staff NGetRestr1 NGetRestr0; % 2012/4/24

INT = 0;

% 健壮性检查
m = prod(q0);
if size(q0,2)~=1
    error('Input variable q must be a s-by-1 vector!\n');
end
if any( q0~=round(q0) ) || any( q0 <= zeros(size(q0)) )
    error('The component of input variable must be positive integer!\n');
end
if beq~=round(beq) || beq <= 0 || beq >= m
    error('Input variable beq must be a positive integer and smaller than prod(q0)!\n');
end
if nargin > 3 + INT && ~isempty(x0) && ( sum(x0==1) ~= beq || sum(x0==1) + sum(x0==0)~= m )
    error('The initial input x0 is wrong!\n');
end


% 全局变量赋初始值 部分1
[H,f] = Fun(q0);
q = q0;
TolXInt = 1e-3;
%options = optimset('display','off','LargeScale','on','Algorithm','interior-point-convex','TolCon',1e-2);
options = optimset('display','off','LargeScale','off','Algorithm','interior-point','TolCon',1e-2);
epsilon = 1e-10;
info1 = zeros(7,1);
n = beq;
% m = length(H);m 之前已赋值
s = length(q);
if nargin > 6 + INT
    opt = opt0;
    clear('opt0');
end
if nargin <= 6 + INT || isempty(opt) || opt(3) < 0 || opt(3) > 1
    opt = zeros(3,1);
    opt(3) = 1;
end
if nargin > 3 +INT && ~isempty(x0) 
    xx = x0;
    clear('x0');
else 
    xx = (randperm(m))';
    v1 = xx(1:beq);
    v0 = xx(beq+1:end);
    xx(v1) = 1;
    xx(v0) = 0;
    clear('v0','v1');
end 

%*************************
xx(:) = 0;
%*************************

% 输入变量不同情况处理
if nargin < 5 + INT || isempty(upperbnd)
    upperbnd = inf;
end
info1(4) = upperbnd;
if nargin < 6 + INT || isempty(lowerbnd)
    lowerbnd = -inf;
end
info1(6) = 1; % 初始化为非 0, 0 表示已达到下界
if lowerbnd <= 0
    info1(6) = -1;
end


% 全局变量赋初始值 部分2
% 此处开始的程序添加于2012/4/12 修改于2012/4/24 2012/4/25
TriesPerLel = zeros(sum(q),1);
NoTriesPerLel = zeros(sum(q),1);
Staff = cumsum([0;q(1:s-1)]);
URestrict = zeros(s,2);
URestrict(:,1) = floor(beq./q);
URestrict(:,2) = mod(beq,q);
URestrict(:,3) = floor((m-beq)./q);
URestrict(:,4) = mod(m-beq,q);%等于q-URestrict(:,2)
NGetRestr1 = zeros(s,1);
NGetRestr0 = zeros(s,1);
% 此处之上程序添加于2012/4/12 修改于2012/4/24 2012/4/25




% a recursive function that processes the BB tree 
% （利用递归法进行二叉树的遍历，实现分枝定界法对整数规划的求解)
bgn = 1; % bgn 表示 rec_Branchbound() 开始寻找可以取 1 的位置
LelCom = ones(s,1) + Staff; % 记录开始寻找可以取 1 的位置 bgn 对应的水平组合
upperbnd = rec_Branchbound(bgn,LelCom,f,beq,upperbnd,lowerbnd);

% 输出汇总
if info1(6) == 0
    % 达到了下界
    info1(5) = lowerbnd;
else    
    info1(5) = upperbnd;
    if info1(6) ~= -1
        info1(6) = (upperbnd-lowerbnd)/lowerbnd;
    end
end
x = xx;
info = info1;
wd2 = -(4/3)^s + 2*info(5)/n^2;
end


%****************************************************************
% a recursive function that processes the BB tree 
% （利用递归法进行二叉树的遍历，实现分枝定界法对整数规划的求解)
%****************************************************************
function upperbnd = rec_Branchbound(bgn,LelCom,f,beq,upperbnd,lowerbnd) 
% 此处程序添加于2012/4/12
global TriesPerLel NoTriesPerLel NGetRestr1 NGetRestr0; 
% URestrict 也是全局变量，JudgeByU1()和JudgeByU0()中用上
% 至此

global H xx epsilon info1 n m s q opt;
global TolXInt options;
info1(3) = info1(3)+1;


if opt(2) %此时将 wd2Ls 中的 neighborhood 技术加入到剪枝技术中
    % 问题规模已降到 4 , 或者 1或 0个数为1 时，采用此块程序   
    if m-bgn+1 <= 4 || min(beq,m-bgn+1-beq) ==1
        v1 = (bgn:bgn+beq-1)'; v0 = (bgn+beq:m)';
        
        
        %**% 2012/5/1 加入此块，解决当 m-bgn+1 == 4，并且 1 的个数为 2
        if m-bgn+1 == 4 && beq ==2
            delta = H(v0(1),v0(2)) + sum(f(v0-bgn+1))-...
                H(v1(1),v1(2))-sum(f(v1-bgn+1));
            if delta < 0
                v = v1;
                v1 = v0;
                v0 = v;
                clear('v');
            end
        end        
        %**% 2012/5/1 块结束
        
        
        delta = 0; % F(son) - F(father)
        for k = 1:beq % 此时 length(v1) == beq
            for l = 1:(m-bgn+1-beq) % 此时 length(v0) == m-bgn+1-beq    
                if length(v1) ~= 1    
                    temp_delta = -( sum(H(v1(k),v1)) - 0.5*H(v1(k),v1(k)) ) + ...
                            ( sum(H(v0(l),[v1(1:k-1);v1(k+1:end);v0(l)])) - ...
                            0.5*H(v0(l),v0(l)) );        
                else
                    temp_delta = 0.5*( H(v0(l),v0(l)) - H(v1,v1) );    
                end
                temp_delta = temp_delta - f(v1(k)-bgn+1) + f(v0(l)-bgn+1);
                if temp_delta < delta    
                    delta = temp_delta;    
                    k_ = k;    
                    l_ = l;    
                end    
            end    
        end
        %if son is better then update i, j, v1, v0    
        if delta < -epsilon              
            temp = v1(k_);    
            v1(k_) = v0(l_);    
            v0(l_) = temp;    
        end
        fval = 0.5*sum(sum(H(v1,v1)));
        if ~isempty(f)
            fval = fval + sum(f(v1-bgn+1));
        end
        if fval < upperbnd - epsilon
            % 此处是第一处可能更新 upperbnd 的位置
            upperbnd = fval;
            xx(v1) = 1;
            xx(v0) = 0;
            info1(1) = info1(1)+1;
        end
        
        if upperbnd <= epsilon + lowerbnd
            % 达到了下界
            info1(6) = 0;
        end
        info1(7) = info1(7) + 1;
        return;
    end    
elseif beq < 1e-8 % beq = 0,说明已经走到叶节点
    if upperbnd > epsilon
        % 此处也是第一处可能更新 upperbnd 的位置
        upperbnd = 0;    
        info1(1) = info1(1)+1;    
    end
    if upperbnd <= epsilon + lowerbnd
        % 达到了下界
        info1(6) = 0;
    end
    info1(7) = info1(7) + 1;
    return;
end
            


if opt(1) % 此时将凸优化加入剪枝手段
    vec = (bgn:m);
    % solve the corresponding QIP model with the integarily constraints removed
    [x,fval,exitflag]=quadprog(H(vec,vec),f,[],[],ones(1,length(vec)),beq,zeros(length(vec),1),ones(length(vec),1),[],options);

    % if the solution is not feasible or the value of the objective function is
    % higher than the current upperbnd return with the input intial solution
    if exitflag<=0 || fval >= upperbnd + epsilon
        return;
    end

    % if the integer-constraint variables turned to be integers within the
    % input tolerance , return
    i = find( abs(x-round(x)) > TolXInt , 1 ); 
    if isempty(i)
        if fval < upperbnd - epsilon    % this solution is better than the current solution hence replace 
        % 此处是第二次可能更新 upperbnd 的位置
            xx(vec) = round(x);     
            upperbnd = fval;
            info1(2) = info1(2)+1;
        end
        if upperbnd <= epsilon + lowerbnd
            % 达到了下界
            info1(6) = 0;
        end
        info1(7) = info1(7) + 1;
        return
    end
    clear('x','fval','vec','i','exitflag');
end



% 2012/4/25 此块程序向前寻找第一个能够满足 U-type 设计的位置
% 向前寻找过程中不断更新 NoTriesPerLel 和 NGetRestr0
% 当从块中跳出时 splt 为寻到的第一个能够满足 U-type 设计的位置
% LelCom 为 splt 对应的水平组合
% splt+1 则应为两个子问题开始寻找第一个能够取1的位置。
NGR0 = zeros(s,1); %用于 NGetRestr0 的回退计算
NTPL = zeros(sum(q),1); % NoTriesPerLelCom 的回退计算
for splt = bgn:m
    [Judge1,Rule1] = JudgeByU1(LelCom);
    if Judge1 
        break; % 找到可以取 1 的位置，跳出 
    end
    % 运行至此处，说明当前 splt 处不能取 1
    [Judge0,Rule0] = JudgeByU0(LelCom);
    if ~Judge0
        % 至此处，说明当前 splt 处也不能取 0, 无解。
        % 这种情形会出现，再继续走下去得到的必然不是 U-type 的。
        % 此时应该将回退 NGetRestr0 和 NoTriesPerLelCom , 再返回。
        if splt > bgn
            NoTriesPerLel = NoTriesPerLel-NTPL;
            NGetRestr0 = NGetRestr0 - NGR0;
        end
        return; 
    end
    NoTriesPerLel(LelCom) = NoTriesPerLel(LelCom) + 1;
    NTPL(LelCom) = NTPL(LelCom)+1;
    NGetRestr0(Rule0) = NGetRestr0(Rule0) + 1;
    NGR0(Rule0) = NGR0(Rule0)+1;
    LelCom = updtLc(LelCom);
end
% 2012/4/25块结束



% 2012/4/25 此块用于决定是否执行取 1 分支子问题
t = min(abs((n-beq+1)/splt-n/m),abs((beq-1)/(m-splt)-n/m));
if t < opt(3)
    C = 0.5*H(splt,splt)+f(splt-bgn+1); % C为进入 1 子问题中需要减去的常数
    % splt 处取 1 分支，则更新NGetRestr1, TriesPerLel。后续要回退。
    NGetRestr1(Rule1) = NGetRestr1(Rule1) + 1;
    TriesPerLel(LelCom) = TriesPerLel(LelCom) + 1;

    % 注意子问题从 splt+1 开始寻找，相应的 LelCom 也应为 splt+1 对应的 LelCom.
    upperbnd1 = rec_Branchbound( splt+1, updtLc(LelCom), f(splt-bgn+2:end,1)+H(splt+1:end,splt), beq-1, upperbnd-C, lowerbnd-C );
    if info1(6) == 0
        % 若取到下界，直接更新 xx 并退出，不需要在更新 upperbnd，
        % 最终的 upperbnd 就是 lowerbnd.
        xx(splt) = 1;
        xx(bgn:splt-1) = 0;    
        return;
    end
    if upperbnd1 < upperbnd - C % if the solution was successfull and gives a better upperbnd
        xx(splt) = 1;
        xx(bgn:splt-1) = 0;
        upperbnd = upperbnd1 + C;
    end
    
    
    %回退 NGetRestr1, TriesPerLel
    NGetRestr1(Rule1) = NGetRestr1(Rule1) - 1;
    TriesPerLel(LelCom) = TriesPerLel(LelCom) - 1;
end
clear('Rule1');
% 2012/4/25 块结束





% 2012/4/25 此块用于决定是否执行 0 分支子问题
% 取 0 分支之前要再次进行判断。判断 splt 处是否可以取 0 分支
[Judge0,Rule0] = JudgeByU0(LelCom);
if ~Judge0
    % 回退 NGetRestr0, NoTriesPerLel
    if splt > bgn
        NoTriesPerLel = NoTriesPerLel-NTPL;
        NGetRestr0 = NGetRestr0 - NGR0;
    end
    return;
end

t = min(abs((n-beq)/splt-n/m),abs(beq/(m-splt)-n/m));
if t < opt(3) 
    % splt 处取 0 分支，则更新NGetRestr0, NoTriesPerLel。后续要回退。
    NoTriesPerLel(LelCom) = NoTriesPerLel(LelCom) + 1;
    NTPL(LelCom) = NTPL(LelCom)+1;
    NGetRestr0(Rule0) = NGetRestr0(Rule0) + 1;
    NGR0(Rule0) = NGR0(Rule0)+1;
    % 注意子问题从 splt+1开始寻找，相应的 LelCom 也应为 splt+1 对应的 LelCom.
    upperbnd0 = rec_Branchbound( splt+1, updtLc(LelCom), f(splt-bgn+2:end), beq, upperbnd, lowerbnd );
    if info1(6) == 0
        % 若取到下界，直接更新 xx 并退出，不需要在更新 upperbnd，
        % 最终的 upperbnd 就是 lowerbnd.
        xx(bgn:splt) = 0;      
        return;    
    end
    if upperbnd0 < upperbnd % if the solution was successfull and gives a better upperbnd
        xx(bgn:splt) = 0;
        upperbnd = upperbnd0;
    end
    NoTriesPerLel = NoTriesPerLel-NTPL;
    NGetRestr0 = NGetRestr0 - NGR0;
end
    
end




% ***************************************************************
% 此函数添加于2012/4/12
% 判断在当前分裂处取 1 是否满足U型设计 Judge1=1表示可以取1
% Rule1记录在分裂处对应因子水平组合中达到均匀的那些水平对应的因子，
% 若分裂取1，则NGetRestr(Rule1) = NGetRestr1(Rule1)+1;
%****************************************************************
function [Judge1,Rule1] = JudgeByU1( LelCom )
global TriesPerLel URestrict NGetRestr1;
if any( TriesPerLel(LelCom) > URestrict(:,1) );
    Judge1 = false;
    Rule1 = [];
    return;
end
Rule1 = TriesPerLel(LelCom) == URestrict(:,1);
% Rule2 = NGetRestr1 >= URestrict(:,2);
% Rule1和Rule2中元素相等并且等于1时才能判断出 Judge1 = 0
Rule = Rule1 & NGetRestr1 >= URestrict(:,2);
if any( Rule )
    Judge1 = false;
    Rule1 = [];
    return;
end
Judge1 = true;
end

% 判断在当前分裂处取 0 是否满足U型设计 Judge0=1表示可以取0
% Rule0记录在分裂处对应因子水平组合中达到均匀的那些水平对应的因子，
% 若分裂取0，则NGetRestr0(Rule1) = NGetRestr1(Rule1)+1;
%****************************************************************
function [Judge0,Rule0] = JudgeByU0( LelCom )
global NoTriesPerLel URestrict NGetRestr0;
if any( NoTriesPerLel(LelCom) > URestrict(:,3) );
    Judge0 = false;
    Rule0 = [];
    return;
end
Rule0 = NoTriesPerLel(LelCom) == URestrict(:,3);
% Rule2 = NGetRestr1 >= URestrict(:,4);
% Rule1和Rule2中元素相等并且等于1时才能判断出 Judge1 = 0
Rule = Rule0 & NGetRestr0 >= URestrict(:,4);
if any( Rule )
    Judge0 = false;
    Rule0 = [];
    return;
end
Judge0 = true;
end


% 2012/4/24 该程序计算LelCom的下一个LelCom，即update LelCom
function LelCom = updtLc(LelCom)
global q s Staff;
LelCom(s) = LelCom(s)+1;
for t = s-1:-1:1
    if LelCom(t+1)-Staff(t+1)>q(t+1)
        LelCom(t) = LelCom(t)+1;
        LelCom(t+1) = LelCom(t+1)-q(t+1);
    else
        break;
    end
end
end

%{
%2012/4/24 该程序计算LelCom的前一个LelCom，即 backward LelCom
function LelCom = backLc(LelCom)
global q s Staff;
LelCom(s) = LelCom(s)-1;
for t = s-1:-1:1
    if LelCom(t+1)-Staff(t+1) < 1
        LelCom(t+1) = LelCom(t+1)+q(t+1);
        LelCom(t) = LelCom(t)-1;
    else
        break;
    end
end
end

% 2012/4/24 该程序计算j位置对应的因子水平组合。
function LelCom = calcLc(bgn)
global q Staff;
% 此程序块添加于2012/4/12
% LelCom：记录即将分裂位置对应的因子水平组合,local,NonPar
% 此块计算要分裂试验对应的水平组合
LelCom = zeros(s,1);
LelCom(s) = mod(bgn,q(s));
if LelCom(s) == 0
    LelCom(s) = q(s);
end
temp = bgn;
for k = s-1:-1:1
    %若按前t个因子的水平组合划分等价类，temp 记录了该等价类的试验序数
    temp = (temp-LelCom(k+1))/q(k+1)+1;
    LelCom(k) = mod(temp,q(k));
    if LelCom(k)==0
        LelCom(k) = q(k);
    end
end
LelCom = LelCom+Staff;
% 块结束
end
%}



