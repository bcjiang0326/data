function [Disc0, Disc_vec, L0] = TA_MD_LHD_1(groups,D,q,OutIter,InIter,T0,T1,Reps,writedata)
% 20150126 
% ���� TA �㷨������ LHD
% ��������ƣ��������� columnwise-pairwise ����
% ������Tang, Xu and Lin (2012) AOS �����uniform fractional factorial design (UFFD). ���ڴ�
%       Tang and Xu (2012) JSPI �����һ�ֹ��� uniform FFD �������㷨�� �����������乹���
%       UFFD ���� OA-based LHD. ������̲��� TA �㷨��
% ׼��Zhou, Fang and Ning (2013) Mixture discrepancy
% ʹ�÷�Χ�����ݸ�����D������ D-based LHD�������Ҫ�� D �� balanced�� 
% Neighborhood: �� D, ( N-by-s design) Ϊ UFFD��D �ĵ�һ��Ϊ q1 ˮƽ���� k = n/q1���������±任
%       0 ---> {0,...,k-1}�� 1 ---> {k,...,2k-1}, �������ơ��������о���˱任�õ�һ�� LHD ��Ϊ L��
%       ��ʱ L �� Neighborhood ����Ϊ �û� L ������ D ��ͬһԪ�ض�Ӧ���еõ������� LHD �ļ���
%       .
% INPUT:
%       D: N-by-n balanced design
%       q: n-by-1 vector, D �����ӵ�ˮƽ��
%       OutIter: �ⲿ��������
%       InIter: �ڲ���������
%       T0: ��ʼ��ֵ��λ��(0,1)֮�䣬����1e-2
%       T1: ������ֵ��λ��(0,1)֮�䣬����1e-6
%       Reps: �ظ����� (default,1)
%       writedata: �Ƿ񽫹�������д���ļ��� 1,д��; 0(default) ��д
% OUTPUT:
%       Disc0: �õ����ž���� Squared Discrepancy
%       Disc_vec: Reps-by-1 vector, ��¼ÿ���ظ������Ž�
%       L0: �õ������ž���

[N,n] = size(D);
if size(q,2)~=1 || n~=size(q,1)
    error(' Input variable q must be a s-by-1 vector!\n');
end
if any(N./q~=floor(N./q))
    error('N and q does not match each other !\n');
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

% �����ļ����
if writedata
    fid = fopen('T_Disc_Disc0.txt','w');
end

% ���ù������
epsilon = 1e-12;

% ���ȼ�����ֵ�仯�ı���,��ÿ InIter �ε�����
% ����һ����ֵ T = T*ratio     
ratio = (T1/T0)^(1/(ceil(OutIter/InIter)-1));

% �����ʼ�� UFFD-based LHD�����ȣ����μ�¼���е� group �������
% groups ��һ��ǰ n/q1 ��Ԫ�ض�Ӧ L ��һ�еĵ�һ�� group���������ơ�
KK = N./q; % ��¼�����У�ÿ��group Ԫ�ظ���
L = zeros(size(D));

Disc0 = inf; 
Disc_vec = zeros(Reps,1);

for rep = 1:Reps
    % �������һ�� L
    for k = 1:n
        for i = 1:q(k)
            v = randperm(KK(k));
            L( groups((i-1)*KK(k)+v,k), k ) = (0:KK(k)-1)'+(i-1)*KK(k);
        end
    end
    
    % ���ȼ����ʼ�� z,sigma1,sigma
    z = (L+0.5)/N-0.5;
    sigma1 = zeros(N,1);
    sigma = zeros(N,N);
    for i = 1:N-1
        for j = i+1:N
            sigma(i,j) = 1/N^2;
            for k = 1:n
                sigma(i,j)=sigma(i,j)*(1.875-0.25*abs(z(i,k))-0.25*abs(z(j,k))-0.75*abs(z(i,k)-z(j,k))+0.5*(z(i,k)-z(j,k))^2);
            end
            sigma(j,i)=sigma(i,j);
        end 
    end
    for i = 1:N
        sigma1(i) = 2/N;
        for k = 1:n
            sigma1(i) = sigma1(i)*(5/3-0.25*abs(z(i,k))-0.25*z(i,k)^2);
        end
        sigma(i,i) = prod(1.875-0.5*abs(z(i,:)))/N^2;
    end
    
    % ��ʼ��Disc
    Disc=0;
    for i = 1:N-1
        for j = i+1:N
            Disc = Disc + sigma(i,j);
        end
    end
    Disc=(19/12)^n+2*Disc+sum(diag(sigma))-sum(sigma1);

    % �жϵ����ظ��ĳ�ʼ���Ƿ����
    Disc_vec(rep) = Disc; 
    if Disc_vec(rep) < Disc0-epsilon
        Disc0 = Disc_vec(rep);
        if nargout > 2
            L0=L;
        end
    end
    
    %�����㷨
    T = T0;
    for Oiter = 1:OutIter
        if mod(Oiter,InIter)==0
            T = T*ratio;
        end
        %����Ҫ�������м���������Ԫ��    
        k=randi(n); % ���ѡ�� k ��
        i0 = randi(q(k)); %���ѡ�����еĵ� i0 ��group
        i1=randi( KK(k) ); % ���ѡ�� group �е� i1 �� i2 ��Ԫ��
        i2=randi( KK(k) );
        while i2==i1    
            i2=randi( KK(k) );
        end
        i = groups( (i0-1)*KK(k)+i1, k);
        j = groups( (i0-1)*KK(k)+i2, k);
        if D(i,k)~=D(j,k)
            error('shit!\n');
        end
    
        %�������ڼ��� delta �� alpha1,alpha2,beta
        alpha1=(5/3-0.25*abs(z(j,k))-0.25*z(j,k)^2)/(5/3-0.25*abs(z(i,k))-0.25*z(i,k)^2);
        alpha2=(1.875-0.5*abs(z(j,k)))/(1.875-0.5*abs(z(i,k)));
        beta = (1.875-0.25*abs(z(j,k))-0.25*abs(z(:,k))-0.75*abs(z(j,k)-z(:,k))+0.5*(z(j,k)-z(:,k)).^2)...
            ./(1.875-0.25*abs(z(i,k))-0.25*abs(z(:,k))-0.75*abs(z(i,k)-z(:,k))+0.5*(z(i,k)-z(:,k)).^2);
    
        %����delta          
        delta = 0;
        for t = 1:N
            if t ~= i && t ~= j
                delta = delta ...
                    +(beta(t)-1)*(sigma(i,t)-sigma(j,t)/beta(t));
            end
        end
        delta = 2*delta+(alpha1-1)*(sigma1(j)/alpha1-sigma1(i))...
            +(alpha2-1)*(sigma(i,i)-sigma(j,j)/alpha2);
        
        %�ж��Ƿ����
        if delta <  Disc*T-epsilon
        %if  T < 1e-5 || delta <  Disc*T-epsilon
            % ���� L,Disc
            temp=L(i,k);L(i,k)=L(j,k);L(j,k)=temp;
            Disc = Disc + delta;
            % ���� z,sigma1,diag(sigma),sigma
            temp=z(i,k);z(i,k)=z(j,k);z(j,k)=temp;            
            sigma1(i)=sigma1(i)*alpha1;
            sigma1(j)=sigma1(j)/alpha1;
            sigma(i,i)=sigma(i,i)*alpha2;
            sigma(j,j)=sigma(j,j)/alpha2;
            for t=1:N
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
                            L0 = L;
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
end

%{
% ���Ժ���
N = 27; n = 5; q = 3*ones(n,1);
X = combins(3,3); % ������
C = [1 1 1; 1 2 0; 1 1 2; 1 0 1; 0 1 2; 1 2 2; 1 1 0]';
B = { [0,1], [0,1,2], [1,1,0,2], [1,1,0,2,1], [1,1,0,2,1,2], [1,1,0,2,1,2,2] };
D0 = mod( [X,X*C(:,1:n-3)+ones(N,1)*B{n-4}], 3);
InIter = 100;
OutIter = 100*InIter;
Reps = 3;
[Disc0, Disc_vec, L0] = TA_MD_LHD(D0,q,OutIter,InIter,1e-3,1e-6,1,0);
fprintf('%.8f\n',Disc0);
fprintf('%.8f\n',MD2_value(L0,ones(n,1)*N));
%}



