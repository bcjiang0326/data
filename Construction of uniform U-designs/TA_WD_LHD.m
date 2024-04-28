function [Disc0, Disc_vec, L0] = TA_WD_LHD(D,q,OutIter,InIter,T0,T1,Reps,writedata,flag)
% 20131001 
% ���� TA �㷨������ LHD
% ��������ƣ��������� columnwise-pairwise ����
% ������Tang, Xu and Lin (2012) AOS �����uniform fractional factorial design (UFFD). ���ڴ�
%       Tang and Xu (2012) JSPI �����һ�ֹ��� uniform FFD �������㷨�� �����������乹���
%       UFFD ���� OA-based LHD. ������̲��� TA �㷨��
% ʹ�÷�Χ�����ݸ�����D������ D-based LHD�������Ҫ�� D �� balanced�� 
% Neighborhood: �� D, ( N-by-n design) Ϊ UFFD��D �ĵ�һ��Ϊ q1 ˮƽ���� k = n/q1���������±任
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
%       flag: 0 ��������(default), 1�������
% OUTPUT:
%       Disc0: �õ����ž���� Squared Discrepancy
%       Disc_vec: Reps-by-1 vector, ��¼ÿ���ظ������Ž�
%       L0: �õ������ž���

[N,n] = size(D);
if size(q,2)~=1 || n~=size(q,1)
    error(' Input variable q must be a n-by-1 vector!\n');
end
if any(N./q~=floor(N./q))
    error('N and q does not match each other !\n');
end
if nargin < 7 || isempty(Reps)
    Reps = 1;
end
if nargin < 8 || isempty(writedata)
    writedata = 0;
end
if nargin < 9 || isempty(flag)
    flag = 0;    
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
groups = zeros(size(D)); 
for j = 1:n
    [~,groups(:,j)] = sort(D(:,j));
    %L(groups(:,j),j) = (1:n)';
end
KK = N./q; % ��¼�����У�ÿ��group Ԫ�ظ���

Disc0 = inf; 
Disc_vec = zeros(Reps,1);

for rep = 1:Reps
    % �������һ�� L
    L = rand_oa_lhd(D,q,flag);
    
    % ���ȼ����ʼ�� z,sigma
    z = L/N;
    sigma = zeros(N,N);
    for i = 1:N-1
        for j = i+1:N
            sigma(i,j) = 1/N^2;
            for k = 1:n
                temp = abs(z(i,k)-z(j,k));
                sigma(i,j)=sigma(i,j)*(1.5-temp*(1-temp));
            end
            sigma(j,i)=sigma(i,j);
        end 
    end
    
    % ��ʼ��Disc
    Disc=0;
    for i = 1:N-1
        for j = i+1:N
            Disc = Disc + sigma(i,j);
        end
    end
    Disc=-(4/3)^n + 1.5^n/N + 2*Disc;

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
    
        %�������ڼ��� delta �� beta
        beta = (1.5-abs(z(j,k)-z(:,k)).*(1-abs(z(j,k)-z(:,k))))...
            ./(1.5-abs(z(i,k)-z(:,k)).*(1-abs(z(i,k)-z(:,k))));
    
        %����delta          
        delta = 0;
        for t = 1:N
            if t ~= i && t ~= j
                delta = delta ...
                    +(beta(t)-1)*(sigma(i,t)-sigma(j,t)/beta(t));
            end
        end
        delta = 2*delta;
        
        %�ж��Ƿ����
        if delta <  Disc*T-epsilon
        %if  (rand()<0.001 && T < 1e-5) || delta <  Disc*T-epsilon
            % ���� L,Disc
            temp=L(i,k);L(i,k)=L(j,k);L(j,k)=temp;
            Disc = Disc + delta;
            % ���� z,sigma1,diag(sigma),sigma
            temp=z(i,k);z(i,k)=z(j,k);z(j,k)=temp;            
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
InIter = 200;
OutIter = 500*InIter;
Reps = 1;
[Disc0, Disc_vec, L0] = TA_WD_LHD(D0,q,OutIter,InIter,1e-3,1e-6,1,0);
fprintf('%.8f\n',Disc0);
%}


