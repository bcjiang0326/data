function [Disc0, Disc_vec, ID0, D0] = TA_REM_WD_f(n,q,OutIter,InIter,T0,T1,Reps,writedata)
% 20130705 
% ������ֵ����(TA)�������ظ��㷨(Repetition Elimination Method)
% The threshold is lowered gemetrically after every InIter itertations
% ��������ƣ��������� columnwise-pairwise ����
% INPUT:
%       n: �������
%       q: s-by-1 vector, ��Ƹ����ӵ�ˮƽ��
%       OutIter: �ⲿ��������
%       InIter: �ڲ���������
%       T0: ��ʼ��ֵ��λ��(0,1)֮�䣬����1e-2
%       T1: ������ֵ��λ��(0,1)֮�䣬����1e-5
%       Reps: �ظ����� (default,1)
%       writedata: �Ƿ񽫹�������д���ļ� -- 1,д��; 0(default)
% OUTPUT:
%       Disc0: �õ����ž���� Squared Discrepancy
%       Disc_vec: Reps-by-1 vector, ��¼ÿ���ظ������Ž�
%       ID0: the rank vector of design D0
%       D0: �õ������ž���
s = size(q,2);
if s~=1
    error('TA_REM_WD: Input variable q must be a s-by-1 vector!\n');
end
if nargin < 7
    Reps = 1;
    writedata = 0;
elseif nargin < 8
    writedata = 0;
end
if Reps > 1
    writedata = 0;
end
% ���ȼ�����ֵ�仯�ı���,��ÿ InIter �ε�����
% ����һ����ֵ T = T*ratio     
ratio = (T1/T0)^(1/(ceil(OutIter/InIter)-1));

s = length(q);
dd = ones(s,1);
for i=1:s-1
    dd(i)=prod(q(i+1:end));
end
N = dd(1)*q(1);
epsilon = 1e-11;

% �������
if writedata
    fid = fopen('T_Disc_Disc0.txt','w');
end

Disc0 = inf; ID0 = [];
Disc_vec = zeros(Reps,1);
for rep = 1:Reps
    % ��������һ��������ظ��ľ������ɳ�ʼ�� sigma
    [D,ID] = rand_U_Type_orth_design(q,n);
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
    Disc=-(4/3)^s+2*Disc+sigma(1,1)*n;
    
    % �жϵ����ظ��ĳ�ʼ���Ƿ����
    Disc_vec(rep) = Disc;
    if Disc_vec(rep) < Disc0-epsilon
        Disc0 = Disc_vec(rep);
        if nargout > 2
            ID0=ID;
        end
    end
    
    Y = zeros(N,1); Y(ID+1) = 1;

    %����TA�㷨
    T = T0;
    for Oiter = 1:OutIter
        if mod(Oiter,InIter)==0
            T = T*ratio;
        end
        %����Ҫ�������м���������Ԫ��    
        isrepeat = 0;
        k=randi(s);
        i=randi(n);
        j=randi(n);
        while j==i || D(i,k)==D(j,k)    
            j=randi(n);
        end
        IDi = ID(i) + (D(j,k)-D(i,k))*dd(k); %�˴���¼�´˶�����
        IDj = ID(j) + (D(i,k)-D(j,k))*dd(k); %���������� ID
       % if ~isempty(find(ID==IDi,1)) || ~isempty(find(ID==IDj,1)) 
       if Y(IDi+1)||Y(IDj+1)
            isrepeat = 1;        
        end   
        %��ʱ���ȿ��ǽ����Ƿ������ظ���
        %�������ظ��У������½��н���
        while isrepeat
            isrepeat = 0;
            k=randi(s);
            i=randi(n);
            j=randi(n);
            while j==i || D(i,k)==D(j,k)
                j=randi(n);
            end
            %�����Ƿ���ܳ����ظ�
            IDi = ID(i) + (D(j,k)-D(i,k))*dd(k);
            IDj = ID(j) + (D(i,k)-D(j,k))*dd(k);
            %if ~isempty(find(ID==IDi,1)) || ~isempty(find(ID==IDj,1)) 
            if Y(IDi+1)||Y(IDj+1)
                isrepeat = 1;        
            end
        end
    
        %�������ڸ��� delta �� beta
        tempj = abs(D(j,k)-D(:,k))/q(k);
        tempi = abs(D(i,k)-D(:,k))/q(k);
        beta = ( 1.5-tempj.*(1-tempj) )./( 1.5-tempi.*(1-tempi) ); 
    
        %����delta          
        delta = 0;
        for t = 1:n
            if t ~= i && t ~= j
                delta = delta ...
                    +(beta(t)-1)*(sigma(i,t)-sigma(j,t)/beta(t));
            end
        end
        delta = delta*2;
        
        %�ж��Ƿ����
        if delta <  Disc*T-epsilon
            % ���� D,ID,Disc,Y
            temp=D(i,k);D(i,k)=D(j,k);D(j,k)=temp;
            Y(ID(i)+1)=0; Y(ID(j)+1)=0; Y(IDi+1)=1; Y(IDj+1)=1;
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
            % �� delta<0 ˵����ʼ�½�.
            if delta < -epsilon
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
            fprintf(fid,'%.8f %.8f %.8f %.8f %.8f\n',delta/Disc,T,Disc,Disc0,Oiter);
        end
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
