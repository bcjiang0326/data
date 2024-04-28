function [Disc0, Disc_vec, ID0, D0] = SA_REM_WD2(n,q,OutIter,InIter,T0,T1,beta_in,Reps,writedata)
% 20130219 
% û�����ȥ�ظ�����
% ��������ƣ��������� columnwise-pairwise ����
% INPUT:
%       n: �������
%       q: s-by-1 vector, ��Ƹ����ӵ�ˮƽ��
%       OutIter: �ⲿ��������
%       InIter: �ڲ���������
%       T0: ��ʼ���¶�ֵ
%       T1: ���յ��¶�ֵ
%       beta_in: �ڲ����²���
%       Reps: �ظ�����
%       writedata: 1,д����; 0,��д(default).
% OUTPUT:
%       Disc0: �õ����ž���� Squared Discrepancy
%       Disc_vec: Reps-by-1 vector, ��¼ÿ���ظ������Ž�
%       ID0: the rank vector of design D0
%       D0: �õ������ž���
s = size(q,2);
if s~=1
    error('SA_REM_WD: Input variable q must be a s-by-1 vector!\n');
end

if nargin < 8
    Reps = 1;
    writedata = 0;
elseif nargin < 9
    writedata = 0;
end
if Reps > 1
    writedata = 0;
end
% �ⲿ���²�������
beta_out = (T1/T0)^(1/(OutIter-1));

s = length(q);
dd = ones(s,1);
for i=1:s-1
    dd(i)=prod(q(i+1:end));
end
epsilon = 1e-11;

%N = dd(1)*q(1);


% �������
if writedata
    fid = fopen('T_Disc_Disc0.txt','w');
end

Disc0 = inf; ID0 = [];
Disc_vec = zeros(Reps,1);
for rep = 1:Reps
    % ��������һ��������ظ��ľ���
    [D,ID] = rand_U_Type_orth_design(q,n);
    
    % ���ɳ�ʼ�� sigma
    sigma = diag((1.5)^s/n^2*ones(n,1));
    for i = 1:n-1
        for j = i+1:n
            sigma(i,j) = 1/n^2;
            for k = 1:s
                temp = abs(D(i,k)-D(j,k))/q(k);
                sigma(i,j)=sigma(i,j)*(1.5-temp*(1-temp));
            end
            sigma(j,i)=sigma(i,j);
        end 
    end
    

    % ��ʼ��Disc
    Disc=0;
    for i = 1:n-1
        for j = i+1:n
            Disc = Disc + sigma(i,j);
        end
    end
    Disc = -(4/3)^s+2*Disc+sigma(1,1)*n;

    % �жϵ����ظ��ĳ�ʼ���Ƿ����
    Disc_vec(rep) = Disc;
    if Disc_vec(rep) < Disc0-epsilon
        Disc0 = Disc_vec(rep);
        if nargout > 2
            ID0=ID;
        end
    end
    
    Tini = T0;
    %����ģ���˻��㷨
    for Oiter = 1:OutIter
        T=Tini;
        ct=0;
        while ct < InIter
            ct = ct + 1;     
            %����Ҫ�������м���������Ԫ��
            k=randi(s);
            i=randi(n);
            j=randi(n);
            while j==i || D(i,k)==D(j,k)    
                j=randi(n);
            end
            IDi = ID(i) + (D(j,k)-D(i,k))*dd(k); %�˴���¼�´˶�����
            IDj = ID(j) + (D(i,k)-D(j,k))*dd(k); %���������� ID
            
            %�������ڸ��� delta �� beta
            tempj = abs(D(j,k)-D(:,k))/q(k);
            tempi = abs(D(i,k)-D(:,k))/q(k);
            beta = ( 1.5-tempj.*(1-tempj) )./( 1.5-tempi.*(1-tempi) ); 
        
            %����delta          
            delta = 0;
            for t = 1:n
                if t ~= i && t ~= j
                    delta = delta ...
                    + (beta(t)-1)*(sigma(i,t)-sigma(j,t)/beta(t));
                end
            end
            delta = delta*2;
        
            %�ж��Ƿ����
            if delta < -epsilon || rand() < exp(-delta/T)
                % ���� D,ID,Disc
                temp=D(i,k);D(i,k)=D(j,k);D(j,k)=temp;
                ID(i) = IDi; ID(j) = IDj;
                Disc = Disc + delta;
                % ���� sigma
                for t=1:n
                    if t~=i && t~=j
                        sigma(i,t)=sigma(i,t)*beta(t);
                        sigma(j,t)=sigma(j,t)/beta(t);
                        sigma(t,i)=sigma(i,t);
                    sigma(t,j)=sigma(j,t);
                    end
                end
                % �� delta<0 ˵����ʼ�½���Ϊʵ�ֳ���½������� ct=0
                if delta < -epsilon
                    ct = 0;
                    % �������ȱ��θ��Ž⣬���¼֮
                    if Disc < Disc_vec(rep) - epsilon
                        Disc_vec(rep) = Disc;
                        % ������ȫ�ָ��Ž⣬���¼֮
                        if Disc_vec(rep) < Disc0-epsilon
                            Disc0 = Disc_vec(rep);
                            if nargout > 2
                                ID0 = ID;
                            end
                        end
                    end
                end
            end
            if writedata
                fprintf(fid,'%.8f %.2f %.8f %.8f %d \n',T,min(exp(-delta/T),1),Disc,Disc_vec(rep),Oiter);
            end
            T = beta_in*T;
        end
        Tini = Tini*beta_out;
    end    
end

if writedata    
    fclose(fid);
end

if nargout > 2
    ID0 = sort(ID0);
end
if nargout > 3
    D0 = Id2Design(q,ID0);
end

end