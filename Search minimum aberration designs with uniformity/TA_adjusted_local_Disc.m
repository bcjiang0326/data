function [Disc0, Disc_vec, D0] = ...
    TA_adjusted_local_Disc(D0,Jset,OutIter,InIter,aa,cc,flag,RR,Reps,outfile)
%20190814 by Bochuan Jiang: adjusted TA algorithm
%�ο�Fang, Ke and Elsawah 2017 ��Construction of uniform designs via
%    an adjusted threshold accepting algorithm��
%ע�⣺�������������ڷǶԳ����
%INPUT:
%    D0: N-by-n matrix, the historical optimal design
%    Jset: k-by-1 vector, k <= n, Ԫ���������У����D0�пɸı����
%    OutIter: �����㷨�е�I���������������ο�Example1������20
%    InIter: �����㷨�е�J���ڲ�����������ο�Example1������5000
%    aa: �����㷨�е�alpha, a scale or Reps-by-1 vector, �ֱ��Ӧ section2.2
%        �� strategy 1 �� strategy 2, �ο�Example1������0.15
%    cc: �����㷨�е�c���ο�Example2������0.03
%    flag:'CD','WD','MD'
%    RR: RR=range(disc)�����ڼ��� T1=aa*RR,
%    Reps: �ظ����� (default,1)
%    outfile: �ļ����������һ���ظ��Ľ��д���ļ� (default,'');
%OUTPUT:
%    Disc0: �õ���������Ƶ�Discֵ
%    Disc_vec: Reps-by-1 vector, ��¼ÿ���ظ������Ž�
%    D0: �õ������ž���

epsilon = 1e-9;
[N,n] = size(D0);
q = zeros(n,1);
for i = 1:n
    q(i) = length(unique(D0(:,i)));
end
if ~issorted(Jset) || any(Jset~=floor(Jset)) || max(Jset) > n || min(Jset) < 1 ||...
        length(unique(Jset))~=length(Jset) || size(Jset,2) ~= 1 || length(size(Jset))>2
    error('Wrong Jset!\n');
end
if OutIter < 1 || InIter < 1
    error('Wrong OutIter, InIter!\n');
end
if ~(0 < cc && cc < 1)
    error('Wrong cc!\n');
end
if nargin < 7 || isempty(flag)
    flag = 'CD';
end
if ~strcmp(flag,'CD') && ~strcmp(flag,'WD') && ~strcmp(flag,'MD')
    error('Wrong flag!\n');
end
if nargin < 8 || isempty(RR)
    RR = Disc_Range(D0,Jset,InIter,flag);
end
if RR < 0
    error('Wrong RR!\n');
end
if nargin < 9 || isempty(Reps) || Reps < 1
    Reps = 1;
end
if length(size(aa))>2 || size(aa,2)~=1 || any(aa<=0) || any(aa>1) ||...
        ( size(aa,1)>1 && size(aa,1)~=Reps )
    error('Wrong aa!\n');
end
if size(aa,1)==1
    aa = aa*ones(Reps,1);
end
if nargin < 10 || isempty(outfile)
    iswrite = 0;
else
    iswrite = 1; % �������
    fid = fopen(outfile,'w');
end

%���ȼ�����ֵ����
T = zeros(OutIter,Reps);
T(1,:) = RR*aa';
for rep = 1:Reps
    for t = 2:OutIter-1
        if T(t-1,rep) > cc*T(1,rep)
            T(t,rep) = (OutIter-t)/OutIter*T(t-1,rep);
        else
            T(t,rep) = (OutIter-t-1)/(OutIter-t)*T(t-1,rep);
        end
    end
end


%��������� historical optimal design ��sigma0,sigma01,Disc0��z0
if strcmp(flag,'MD')
    [Disc0,sigma0,sigma01] = MD2_value(D0,q);
    z0 = (D0+0.5)./(ones(N,1)*q')-0.5;
elseif strcmp(flag,'CD')
    [Disc0,sigma0,sigma01] = CD2_value(D0,q);
    z0 = 0.5*( (D0+0.5)./(ones(N,1)*q')-0.5 );
elseif strcmp(flag,'WD')
    [Disc0,sigma0] = WD2_value(D0,q);
end

%��¼�����ظ����õ������Ž�
Disc_vec = zeros(Reps,1);

%�����㷨
D = D0; %current design
for rep = 1:Reps
    % ���ȶ�Jset��ǵ��н������������Ϊ��rep�ظ��ĳ�ʼ����
    for j = Jset'
        D(:,j) = D0(randperm(N,N),j);
    end
    % ���ɳ�ʼ��sigma,sigma1,Disc �� z
    if strcmp(flag,'MD')
        [Disc,sigma,sigma1] = MD2_value(D,q);
        z = (D+0.5)./(ones(N,1)*q')-0.5;
    elseif strcmp(flag,'CD')
        [Disc,sigma,sigma1] = CD2_value(D,q);
        z = 0.5*( (D+0.5)./(ones(N,1)*q')-0.5 );
    elseif strcmp(flag,'WD')
        [Disc,sigma] = WD2_value(D,q);
    end
    
    %�����ظ��ĳ�ʼ����ȫ���ظ��и��ţ��ǣ�����D0,sigma0,sigma01
    Disc_vec(rep) = Disc;
    if Disc_vec(rep) < Disc0-epsilon
        D0 = D;
        Disc0 = Disc_vec(rep);
        sigma0 = sigma;
        if strcmp(flag,'CD') || strcmp(flag,'MD')
            sigma01 = sigma1;
            z0 = z;
        end
    end
    
    %�������ѭ��
    for Oiter = 1:OutIter
        %Current Design����D0���ǣ�ȡD0��ΪCurrent Design
        if Disc > Disc0+epsilon
            D = D0;
            Disc = Disc0;
            sigma = sigma0;
            if strcmp(flag,'CD') || strcmp(flag,'MD')
                sigma1 = sigma01;
                z = z0;
            end
            %�µ�Current Design�ڱ����ظ��и��ţ��ǣ����±����ظ����Ž�
            if Disc < Disc_vec(rep) - epsilon
                Disc_vec(rep) = Disc;
            end
        end
        
        %�����ڲ�ѭ��
        for Iiter = 1:InIter
            %����Ҫ�������м���������Ԫ��
            k = Jset(randi(length(Jset)));
            i = randi(N);
            j = randi(N);
            while j==i || abs(D(i,k)-D(j,k))<epsilon
                j=randi(N);
            end
            
            %�������ڸ��� delta �� alpha1,alpha2,beta
            if strcmp(flag,'MD')
                alpha1=(5/3-0.25*abs(z(j,k))-0.25*z(j,k)^2)/(5/3-0.25*abs(z(i,k))-0.25*z(i,k)^2);
                alpha2=(1.875-0.5*abs(z(j,k)))/(1.875-0.5*abs(z(i,k)));
                beta = (1.875-0.25*abs(z(j,k))-0.25*abs(z(:,k))-0.75*abs(z(j,k)-z(:,k))+0.5*(z(j,k)-z(:,k)).^2)...
                    ./(1.875-0.25*abs(z(i,k))-0.25*abs(z(:,k))-0.75*abs(z(i,k)-z(:,k))+0.5*(z(i,k)-z(:,k)).^2);
            elseif strcmp(flag,'CD')
                alpha1=(1+abs(z(j,k))-2*z(j,k)^2)/(1+abs(z(i,k))-2*z(i,k)^2);
                alpha2=(1+2*abs(z(j,k)))/(1+2*abs(z(i,k)));
                beta = (1+abs(z(j,k))+abs(z(:,k))-abs(z(j,k)-z(:,k)))...
                    ./(1+abs(z(i,k))+abs(z(:,k))-abs(z(i,k)-z(:,k)));
            elseif strcmp(flag,'WD')
                tempj = abs(D(j,k)-D(:,k))/q(k);
                tempi = abs(D(i,k)-D(:,k))/q(k);
                beta = ( 1.5-tempj.*(1-tempj) )./( 1.5-tempi.*(1-tempi) );
            end
            
            %����delta
            delta = 0;
            for t = 1:N
                if t ~= i && t ~= j
                    delta = delta ...
                        +(beta(t)-1)*(sigma(i,t)-sigma(j,t)/beta(t));
                end
            end
            delta = 2*delta;
            if strcmp(flag,'CD') || strcmp(flag,'MD')
                delta = delta +(alpha1-1)*(sigma1(j)/alpha1-sigma1(i))...
                    +(alpha2-1)*(sigma(i,i)-sigma(j,j)/alpha2);
            end
            
            %�ж��Ƿ����
            if delta < T(Oiter,rep)-epsilon
                % ���� D,Disc
                temp=D(i,k);D(i,k)=D(j,k);D(j,k)=temp;
                Disc = Disc + delta;
                % ���� sigma �ǶԽǲ���
                for t=1:N
                    if t~=i && t~=j
                        sigma(i,t)=sigma(i,t)*beta(t);
                        sigma(j,t)=sigma(j,t)/beta(t);
                        sigma(t,i)=sigma(i,t);
                        sigma(t,j)=sigma(j,t);
                    end
                end
                if strcmp(flag,'CD') || strcmp(flag,'MD')
                    % ���� z,sigma1,diag(sigma)
                    temp=z(i,k);z(i,k)=z(j,k);z(j,k)=temp;
                    sigma1(i)=sigma1(i)*alpha1;
                    sigma1(j)=sigma1(j)/alpha1;
                    sigma(i,i)=sigma(i,i)*alpha2;
                    sigma(j,j)=sigma(j,j)/alpha2;
                end
                
                % �� delta<0 ˵����ʼ�½�.
                if delta < -epsilon
                    % �������ȱ��θ��Ž⣬���¼֮
                    if Disc < Disc_vec(rep) - epsilon
                        Disc_vec(rep) = Disc;
                        % ������ȫ�ָ��Ž⣬���¼֮
                        if Disc < Disc0-epsilon
                            D0 = D;
                            Disc0 = Disc;
                            sigma0 = sigma;
                            if strcmp(flag,'CD') || strcmp(flag,'MD')
                                sigma01 = sigma1;
                                z0 = z;
                            end
                        end
                    end
                end
            end
            if iswrite && rep == Reps
                fprintf(fid,'%.8f\n',Disc);
            end
        end
    end
end

if iswrite
    fclose(fid);
end
end


%{
%���Դ���
N = 27; n = 13; s = 3;
q = s*ones(n,1);
D0 = rand_U_Type_design(q,N);
Jset = (1:n)';

InIter = 5000;
OutIter = 20;
aa = [0.15,0.016,0.01,0.002,0.0005]';
cc = 0.03;
flag = 'WD';
Reps = 1;
outfile = 'out.txt';

figure
for i = 1:length(aa)
    [Disc1, Disc_vec, D0] = TA_adjusted_local_Disc(D0,Jset,OutIter,InIter,aa(i),cc,flag,'',Reps,outfile);
    if strcmp(flag,'MD')
        Disc = MD2_value(D0);
    elseif strcmp(flag,'CD')
        Disc = CD2_value(D0);
    elseif strcmp(flag,'WD')
        Disc = WD2_value(D0);
    end
    fprintf('%.6f  %.6f  %.6f  %d\n', Disc1, Disc, mean(Disc_vec),isOA(D0,2));
    aveDisc_seq = importdata('out.txt');    
    plot(1:length(aveDisc_seq),aveDisc_seq);
    hold on
end
hold off
%}
