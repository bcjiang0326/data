function [Disc0, Disc_vec, ID0, D0] = SA_REM_CD(n,q,OutIter,InIter,T0,T1,beta_in,Reps,writedata)
% 20121203 
% ����ģ���˻�(SA)�������ظ��㷨(Repetition Elimination Method)
% ��������ƣ��������� columnwise-pairwise ����
% INPUT:
%       n: �������
%       q: s-by-1 vector, ��Ƹ����ӵ�ˮƽ��
%       OutIter: �ⲿ��������
%       InIter: �ڲ���������
%       T0: ��ʼ���¶�
%       T1: ���յ��¶�
%       beta_in: �ڲ����²���
%       Reps: �ظ�����
%       writedata: 1,д����; 0,��д(default).
% OUTPUT:
%       disc: �õ����ž���� Squared Discrepancy
%       Disc_vec: Reps-by-1 vector, ��¼ÿ���ظ������Ž�
%       ID0: the rank vector of design D0
%       D0: �õ������ž���



s = size(q,2);
if s~=1
    error('SA_REM: Input variable q must be a s-by-1 vector!\n');
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

N = dd(1)*q(1);

if writedata
    fid = fopen('T_Disc_Disc0.txt','w');
end

Disc0 = inf; ID0 = [];
Disc_vec = zeros(Reps,1);
for rep = 1:Reps
    % ��������һ��������ظ��ľ���
    [D,ID] = rand_U_Type_orth_design(q,n);
    Y = zeros(N,1); Y(ID+1) = 1;
    %���ɳ�ʼ�� z,sigma1,sigma
    z = zeros(n,s);
    sigma1 = zeros(n,1);
    sigma = zeros(n,n);
    for i = 1:n
        z(i,:) = 0.5*( (D(i,:)+0.5)./q'-0.5 );
    end
    for i = 1:n-1
        for j = i+1:n
            sigma(i,j) = 1/n^2;
            for k = 1:s
                sigma(i,j)=sigma(i,j)*(1+abs(z(i,k))+abs(z(j,k))-abs(z(i,k)-z(j,k)));
            end
            sigma(j,i)=sigma(i,j);
        end 
    end
    for i = 1:n
        sigma1(i) = 2/n;
        for k = 1:s
            sigma1(i) = sigma1(i)*(1+abs(z(i,k))-2*z(i,k)^2);
        end
        sigma(i,i) = prod(1+2*abs(z(i,:)))/n^2;
    end
    % ��ʼ��Disc
    Disc=0;
    for i = 1:n-1
        for j = i+1:n
            Disc = Disc + sigma(i,j);
        end
    end
    Disc=(13/12)^s+2*Disc+sum(diag(sigma))-sum(sigma1);
    % ��ʼ�����Ž���Ϣ
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
            isrepeat = 0;
            k=randi(s);
            i=randi(n);
            j=randi(n);
            while j==i || D(i,k)==D(j,k)    
                j=randi(n);
            end
            IDi = ID(i) + (D(j,k)-D(i,k))*dd(k); %�˴���¼�´˶�����
            IDj = ID(j) + (D(i,k)-D(j,k))*dd(k); %���������� ID
            if Y(IDi+1) || Y(IDj+1)
            %if ~isempty(find(ID==IDi,1)) || ~isempty(find(ID==IDj,1))        
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
                if Y(IDi+1) || Y(IDj+1)
                    isrepeat = 1;
                end
            end
            %�������ڼ��� delta �� alpha1,alpha2,beta, ������ delta
            alpha1=(1+abs(z(j,k))-2*z(j,k)^2)/(1+abs(z(i,k))-2*z(i,k)^2);
            alpha2=(1+2*abs(z(j,k)))/(1+2*abs(z(i,k)));
            beta = (1+abs(z(j,k))+abs(z(:,k))-abs(z(j,k)-z(:,k)))...
                ./(1+abs(z(i,k))+abs(z(:,k))-abs(z(i,k)-z(:,k)));
            delta = 0;
            for t = 1:n
                if t ~= i && t ~= j
                    delta = delta ...
                        + (beta(t)-1)*(sigma(i,t)-sigma(j,t)/beta(t));
                end
            end
            delta = 2*delta+(alpha1-1)*(sigma1(j)/alpha1-sigma1(i))...
                +(alpha2-1)*(sigma(i,i)-sigma(j,j)/alpha2);
        
            %�ж��Ƿ����
            if delta < -epsilon || rand() < exp(-delta/Disc/T)
                % ���� D,ID,Disc
                Y(ID(i)+1) = 0; Y(ID(j)+1) = 0;
                Y(IDi+1) = 1; Y(IDj+1) = 1;
                temp = D(i,k); D(i,k) = D(j,k); D(j,k) = temp;
                ID(i) = IDi; ID(j) = IDj;
                Disc = Disc + delta;
                % ���� z,sigma1,diag(sigma),sigma
                temp=z(i,k);z(i,k)=z(j,k);z(j,k)=temp;            
                sigma1(i)=sigma1(i)*alpha1;
                sigma1(j)=sigma1(j)/alpha1;
                sigma(i,i)=sigma(i,i)*alpha2;
                sigma(j,j)=sigma(j,j)/alpha2;
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
                fprintf(fid,'%.8f %.2f %.8f %.8f %d\n',T,min(exp(-delta/Disc/T),1),Disc,Disc0,Oiter);
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