function [A0,A_vec,D0] = GMA(N,n,s,Iter,Reps,writedata)
% 20131001 
% ���� greedy with mutation �㷨���Գ� GMA ���
% �������� columnwise-pairwise ���ԣ��������̲��ܱ�֤Ϊ OA
% INPUT:
%       N: run-size
%       n: number of factors
%       s: number of levels
%       Iter: �ܵĵ�������
%       Reps: �ظ����� (default,1)
%       writedata: �Ƿ񽫹�������д���ļ��� 1,д��; 0(default) ��д
% OUTPUT:
%       A0: �õ����ž���� GWP
%       A_vec: 
%       D0: �õ������ž���

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

% �����ļ����
if writedata
    fid = fopen('GMA.txt','w');
end

% ���ù������
epsilon = 1e-12;


% �������ɲ�����
q = s*ones(n,1);
P = Krawtchouk_polyMat(n,s);

% �����ʼ GWP ֵ
A0 = zeros(1,n); A0(1) = inf;
A_vec = zeros(Reps,n); A_vec(:,1) = inf;

for rep = 1:Reps
    % ��������һ�����ƽ�����
    D = rand_U_Type_design(q,N);
    %j �����ʼ�� distance distribution
    [B,MatH] = dist_distr(D);
    % �����ʼ�� GWP
    A= B*P/N;
    
    % ����������µĵ��Ͻ�ֵ
    maxT = 0.1*A; 
    
    
     % �жϵ����ظ��ĳ�ʼ���Ƿ����
     A_vec(rep,:) = A;
     v = find(abs(A-A0)>epsilon,1);
     if A(v) < A0(v)
         A0 = A;
         if nargout > 2
             D0 = D;
         end
     end
     
     % �����㷨
     for Oiter = 1:Iter
          %����Ҫ�������м���������Ԫ��    
          k0 = randi( n ); % ���ѡ�� k0 ��
          i0 = randi( N );
          j0 = randi( N );
          while j0==i0 || D(i0,k0) == D(j0,k0) % || MatH(i0,j0)==1 %ʹ���û���� D�������ظ�
              j0=randi( N );
          end
          
          %���� delta_B
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
          
          % ���� delta_A
          v = dt_B>epsilon | dt_B<-epsilon;
          dt_A = dt_B(v)/N*P(v,:);
          
          
          % �ж��Ƿ����
          v = find(abs(dt_A)>epsilon,1);
          %if  isempty(v) || (dt_A(v) < 0) || (rand()<0.01 && dt_A(v) <  A(v)*T-epsilon)
          if  isempty(v) || (dt_A(v) < 0) || ( rand()<0.05 && dt_A(v)<maxT(v)*0.1 )
                % ���� A,B,MatH,D
                temp=D(i0,k0); D(i0,k0)=D(j0,k0); D(j0,k0)=temp;
                MatH(i0,ti) = MatH(i0,ti)+1; MatH(ti,i0) = MatH(i0,ti)';
                MatH(i0,tj) = MatH(i0,tj)-1; MatH(tj,i0) = MatH(i0,tj)';
                MatH(j0,ti) = MatH(j0,ti)-1; MatH(ti,j0) = MatH(j0,ti)';
                MatH(j0,tj) = MatH(j0,tj)+1; MatH(tj,j0) = MatH(j0,tj)';
                B = B+dt_B;
                A = A+dt_A;
                % �� delta_A < 0 ˵����ʼ�½�
                if dt_A(v)<-epsilon
                    % �������ȱ��θ��Ž⣬���¼֮
                    v = find(abs(A-A_vec(rep,:))>epsilon,1);
                    if ~isempty(v) && A(v) < A_vec(rep,v)
                        A_vec(rep,:) = A;
                        % ������ȫ�ָ��Ž⣬���¼֮
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
% ���Ժ���
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

