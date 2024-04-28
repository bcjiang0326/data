function [A0,A_vec,D0] = GMA(N,n,s,Iter,Reps,writedata)
% 20131001 
% 基于 greedy with mutation 算法求解对称 GMA 设计
% 迭代采用 columnwise-pairwise 策略，迭代过程不能保证为 OA
% INPUT:
%       N: run-size
%       n: number of factors
%       s: number of levels
%       Iter: 总的迭代次数
%       Reps: 重复次数 (default,1)
%       writedata: 是否将过程数据写入文件。 1,写入; 0(default) 不写
% OUTPUT:
%       A0: 得到最优矩阵的 GWP
%       A_vec: 
%       D0: 得到的最优矩阵

if  N ~= floor(N) || n ~= floor(n) || s ~= floor(s) || N/s ~= floor(N/s) 
    error('Input errors!\n');
end
if nargin < 5
    Reps = 1;
    writedata = 0;
elseif nargin < 6
    writedata = 0;
end
if Reps > 1
    writedata = 0;
end

% 设置文件输出
if writedata
    fid = fopen('GMA.txt','w');
end

% 设置归零参数
epsilon = 1e-12;


% 计算若干不变量
q = s*ones(n,1);
P = Krawtchouk_polyMat(n,s);

% 给予初始 GWP 值
A0 = zeros(1,n); A0(1) = inf;
A_vec = zeros(Reps,n); A_vec(:,1) = inf;

for rep = 1:Reps
    % 首先生成一个随机平衡矩阵
    D = rand_U_Type_design(q,N);
    %j 计算初始的 distance distribution
    [B,MatH] = dist_distr(D);
    % 计算初始的 GWP
    A= B*P/N;
    
    % 给定允许更新的的上界值
    maxT = 0.1*A; 
    
    
     % 判断当次重复的初始解是否更优
     A_vec(rep,:) = A;
     v = find(abs(A-A0)>epsilon,1);
     if A(v) < A0(v)
         A0 = A;
         if nargout > 2
             D0 = D;
         end
     end
     
     % 进入算法
     for Oiter = 1:Iter
          %产生要交换的列及列中两个元素    
          k0 = randi( n ); % 随机选第 k0 列
          i0 = randi( N );
          j0 = randi( N );
          while j0==i0 || D(i0,k0) == D(j0,k0) % || MatH(i0,j0)==1 %使得置换后的 D不能有重复
              j0=randi( N );
          end
          
          %计算 delta_B
          dt_B  = zeros(1,n+1);
          ti = D(:,k0)==D(i0,k0); ti(i0) = 0;
          tj = D(:,k0)==D(j0,k0); tj(j0) = 0;
          for a = MatH(i0,ti)
              dt_B(a+1) = dt_B(a+1)-2/N;
              dt_B(a+2) = dt_B(a+2)+2/N;
          end
          for a = MatH(j0,tj)
              dt_B(a+1) = dt_B(a+1)-2/N;
              dt_B(a+2) = dt_B(a+2)+2/N;
          end
          for a = MatH(i0,tj)
              dt_B(a+1) = dt_B(a+1)-2/N;
              dt_B(a) = dt_B(a)+2/N;
          end
          for a = MatH(j0,ti)
              dt_B(a+1) = dt_B(a+1)-2/N;
              dt_B(a) = dt_B(a)+2/N;
          end
          
          % 计算 delta_A
          v = dt_B>epsilon | dt_B<-epsilon;
          dt_A = dt_B(v)/N*P(v,:);
          
          
          % 判断是否更新
          v = find(abs(dt_A)>epsilon,1);
          %if  isempty(v) || (dt_A(v) < 0) || (rand()<0.01 && dt_A(v) <  A(v)*T-epsilon)
          if  isempty(v) || (dt_A(v) < 0) || ( rand()<0.05 && dt_A(v)<maxT(v)*0.1 )
                % 更新 A,B,MatH,D
                temp=D(i0,k0); D(i0,k0)=D(j0,k0); D(j0,k0)=temp;
                MatH(i0,ti) = MatH(i0,ti)+1; MatH(ti,i0) = MatH(i0,ti)';
                MatH(i0,tj) = MatH(i0,tj)-1; MatH(tj,i0) = MatH(i0,tj)';
                MatH(j0,ti) = MatH(j0,ti)-1; MatH(ti,j0) = MatH(j0,ti)';
                MatH(j0,tj) = MatH(j0,tj)+1; MatH(tj,j0) = MatH(j0,tj)';
                B = B+dt_B;
                A = A+dt_A;
                % 当 delta_A < 0 说明开始下降
                if dt_A(v)<-epsilon
                    % 若遇到比本次更优解，则记录之
                    v = find(abs(A-A_vec(rep,:))>epsilon,1);
                    if ~isempty(v) && A(v) < A_vec(rep,v)
                        A_vec(rep,:) = A;
                        % 若遇到全局更优解，则记录之
                        v = find(abs(A-A0)>epsilon,1);
                        if (~isempty(v) &&A(v) < A0(v))
                            A0 = A;
                            if nargout > 2
                                D0 = D;
                            end
                        end
                    end
                end
          end
          if writedata
              fprintf(fid,'%.8f ',dt_A(2:end)./maxT(2:end));
              fprintf(fid,'\n');
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
A = GWP(D,s);

InIter = 100;
OutIter = 100*InIter;
Reps = 1;
[A0, A_vec, D0] = GMA(N,n,s,OutIter,Reps,0);
fprintf('%.2f ',A0);fprintf('\n');
AA = GWP(D0,s);
norm(A0-AA);

data = importdata('GMA.txt');
data = unique(data,'rows');
%}

