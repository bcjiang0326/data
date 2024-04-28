function [Disc0, Disc_vec, ID0, D0] = SA_IPM(n,q,OutIter,InIter,T0,beta_in,beta_out,isWD,Reps,writedata)
% 20121230 “ZhouFangNing2012"
% 基于模拟退火整数规划(SA)的消除重复算法(Repetition Elimination Method)
% 求解均匀设计，迭代采用 columnwise-pairwise 策略
% INPUT:
%       n: 试验次数
%       q: s-by-1 vector, 设计各因子的水平数
%       OutIter: 外部迭代次数
%       InIter: 内部迭代次数
%       T0: 初始的温度
%       beta_in: 内部降温参数
%       beta_out: 外部降温参数
%       isWD: 0,CD; 1,WD.
%       Reps: 重复次数 1(default)
%       writedata: 1,写数据; 0,不写(default).
% OUTPUT:
%       disc: 得到最优矩阵的 Squared Discrepancy 
%       D0: 得到的最优矩阵

if min(size(q)) ~= 1
    error('SA_REM_WD: Input variable q must be a s-by-1 vector!\n');
end
if nargin < 8
    isWD = 1;
    Reps = 1;
    writedata = 0;
elseif nargin < 9
    Reps = 1;
    writedata = 0;
elseif nargin < 10
    writedata = 0;
end
if Reps > 1
    writedata = 0;
end

s = length(q);
N = prod(q);
epsilon = 1e-11;

% 产生初始迭代用的矩阵
if isWD
    Q = getMatWD(q);    
else
    [Q,b] = getMatCD(q);
    Q = Q-ones(N,1)*b'-b*ones(1,N);
end
%c = min(min(Q));Q = Q-c;KA = sum(sum(Q))+1;
KA = 2*sum(sum(abs(Q)))+1; 
Q = Q/KA+1-2*n*eye(N);

% 设置输出
if writedata
    fid = fopen('T_Disc_Disc0.txt','w');
end


% 最优解设置
ObjVal0 = inf;
ID0 = [];
Disc_vec = zeros(Reps,1);
for rep = 1:Reps
    %产生初始的解 y, ObjVal, gains g;
    y = rand(N,1)<0.5;
    ObjVal = sum(sum(Q(y,y)));
    Disc_vec(rep) = ObjVal;
    if ObjVal < ObjVal0-epsilon
        ObjVal0 = ObjVal;
        if nargout > 2
            ID0 = find(y);
        end
    end
    % calculate all gains g_i of y
    g = zeros(N,1);
    for i = 1:N
        g(i) = Q(i,i)+(1-2*y(i))*2*sum(Q(y,i));
    end
    
    Tini = T0;
    %进入模拟退火算法外层循环
    for Oiter = 1:OutIter
        T=Tini;
        %进入模拟退火算法内层循环
        ct=0;
        while ct < InIter
            ct = ct + 1;
            %此处要做 N 次内层循环。。。
            for t = 1:N
                [nabla,j] = min(g);
                if nabla < 0 %说明此时找到了一个更好的解
                    ct = 0;
                    ObjVal = ObjVal + g(j);
                    y(j) = ~y(j); % the bit j is flipped
                    %update all gains g_i                    
                    temp = -g(j);
                    g = g+2*Q(:,j).*(1-2*y)*(2*y(j)-1);
                    g(j)=temp;
                    
                    if ObjVal < Disc_vec(rep)-epsilon
                        Disc_vec(rep) = ObjVal;
                        if Disc_vec(rep) < ObjVal0-epsilon
                            ObjVal0 = Disc_vec(rep);
                            if nargout > 2
                                y0 = y;
                            end
                        end
                    end
                elseif rand() < exp(-g(j)/T)
                    j = randi(N); % randomly choose j \in {1,...,N}, flip bit j.
                    y(j) = ~y(j); % the bit j is flipped
                    ObjVal = ObjVal + g(j);
                    %update all gains g_i
                    temp = -g(j);
                    g = g+2*Q(:,j).*(1-2*y)*(2*y(j)-1);
                    g(j)=temp;
                end
                if writedata
                    fprintf(fid,'%.8f %.2f %.8f %.8f %d \n',T,min(exp(-g(j)/T),1),ObjVal,Disc_vec(rep),Oiter);
                end
            end
            T = beta_in*T;
        end
        Tini = Tini*beta_out;       
    end
end

if nargout > 2
    ID0 = find(y0);
    ID0 = ID0-1;
end
if nargout > 3
    D0 = Id2Design(q,ID0);
end

if isWD
    Disc0 = -(4/3)^s + KA*(ObjVal0/n^2+1);
    if nargout > 1
        Disc_vec = -(4/3)^s + KA*(Disc_vec/n^2+1);
    end
else
    Disc0 = (13/12)^s + KA*(ObjVal0/n^2+1);
    if nargout > 1
        Disc_vec = (13/12)^s + KA*(Disc_vec/n^2+1);
    end
end

end

%{
% test SA_IPM
% ZhouFangNing2012 搜索设计SA_IPM
levels = 3; s = 5; n = 48;

q = levels*ones(s,1);
N = prod(q);

% ZhouFangNing2012 输入参数
InIter = 10;
T0 = 1/q(1);
beta_in = 0.99;
Reps = 1;
isWD = 1;
if N < 500
    OutIter = 10; beta_out = 0.9;
else
    OutIter = 2; beta_out = 0.8;
end
t = cputime;
[Disc0, Disc_vec, ID0, D0] = SA_IPM(n,q,OutIter,InIter,T0,beta_in,beta_out,isWD,Reps);
cout_8f(Disc0); cout_8f(cputime-t);
%}
