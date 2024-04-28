function [Disc0, Disc_vec, ID0, D0] = SA_IPM(n,q,OutIter,InIter,T0,beta_in,beta_out,isWD,Reps,writedata)
% 20121230 ��ZhouFangNing2012"
% ����ģ���˻������滮(SA)�������ظ��㷨(Repetition Elimination Method)
% ��������ƣ��������� columnwise-pairwise ����
% INPUT:
%       n: �������
%       q: s-by-1 vector, ��Ƹ����ӵ�ˮƽ��
%       OutIter: �ⲿ��������
%       InIter: �ڲ���������
%       T0: ��ʼ���¶�
%       beta_in: �ڲ����²���
%       beta_out: �ⲿ���²���
%       isWD: 0,CD; 1,WD.
%       Reps: �ظ����� 1(default)
%       writedata: 1,д����; 0,��д(default).
% OUTPUT:
%       disc: �õ����ž���� Squared Discrepancy 
%       D0: �õ������ž���

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

% ������ʼ�����õľ���
if isWD
    Q = getMatWD(q);    
else
    [Q,b] = getMatCD(q);
    Q = Q-ones(N,1)*b'-b*ones(1,N);
end
%c = min(min(Q));Q = Q-c;KA = sum(sum(Q))+1;
KA = 2*sum(sum(abs(Q)))+1; 
Q = Q/KA+1-2*n*eye(N);

% �������
if writedata
    fid = fopen('T_Disc_Disc0.txt','w');
end


% ���Ž�����
ObjVal0 = inf;
ID0 = [];
Disc_vec = zeros(Reps,1);
for rep = 1:Reps
    %������ʼ�Ľ� y, ObjVal, gains g;
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
    %����ģ���˻��㷨���ѭ��
    for Oiter = 1:OutIter
        T=Tini;
        %����ģ���˻��㷨�ڲ�ѭ��
        ct=0;
        while ct < InIter
            ct = ct + 1;
            %�˴�Ҫ�� N ���ڲ�ѭ��������
            for t = 1:N
                [nabla,j] = min(g);
                if nabla < 0 %˵����ʱ�ҵ���һ�����õĽ�
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
% ZhouFangNing2012 �������SA_IPM
levels = 3; s = 5; n = 48;

q = levels*ones(s,1);
N = prod(q);

% ZhouFangNing2012 �������
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
