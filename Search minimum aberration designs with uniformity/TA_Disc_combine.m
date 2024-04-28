function [Disc0, Disc0vec, D0] = TA_Disc_combine(D,groups,AllIter,InIter,T0,T1,flag,Reps,filename)
%20190822 by Bochuan Jiang
%�ο� Jiang and Ai (2019) New construction of uniform designs via average
%     discrepancy theory
%����ư��зֳ���������ƣ��ڵ������������ѡ��һ������ƣ�����û��������У�����ʵ��ƫ����С����
%INPUT:
%   D: ������ơ�
%   groups: m-by-1 cell, ��ư��з�Ϊm�飬����group{i} Ϊnk-by-1 vector, ��ʾ��i���б�š�
%   AllIter: �ܵ������� = �ڲ��������� * �ⲿ����������
%   InIter: �ڲ�����������
%   T0: ��ʼ��ֵ��λ��(0,1)֮�䣬����1e-2��
%   T1: ������ֵ��λ��(0,1)֮�䣬����1e-5��
%   flag: 'CD','WD' or 'MD'��
%   Reps: �ظ��������� (default,1)��
%   filename: �ṩ�ļ����Ļ����򽫽��д���ļ���
%OUTPUT:
%   Disc0: ����ƫ��ֵ
%   Disc0vec: Reps�ظ������õ�������ƫ��ֵ����
%   D0: MAD design

epsilon = 1e-10;
[N,n] = size(D);
q = zeros(n,1);
for i = 1:n
    q(i) = length(unique(D(:,i)));
end
g = max(size(groups));
if length(size(groups,2))> 2 || min(size(groups))~=1
    error('Wrong groups!\n');
end
q_ = zeros(0,1);
for k = 1:g
    q_ = cat(1,q_,groups{k});
end
if any((1:n)'~=sort(q_))
    error('Wrong groups!\n');
end
clear('q_');

if AllIter < 1 || InIter < 1 || AllIter < InIter
    error('TA_MAD: Wrong OutIter, InIter!\n');
end
if ~(0 <= T1 && T1 <= T0 && T0 < 1)
    error('TA_MAD: Wrong T0, T1!\n');
end
if ~strcmp(flag,'CD') && ~strcmp(flag,'WD') && ~strcmp(flag,'MD')
    error('Wrong flag!\n');
end
if nargin < 8 || isempty(Reps) || Reps < 1
    Reps = 1;
end
if nargin < 9 || isempty(filename)
    iswrite = 0;
else
    iswrite = 1; % �������
    fid = fopen(filename,'w');
    rep0 = 1; %��¼��rep0����������
end

%���ȼ������������ D �йز���
if strcmp(flag,'MD') || strcmp(flag,'CD')
    [Disc,sigma,sigma1] = Disc2_value(D,q,flag);
    sgm_k1 = zeros(N,g);
    for k = 1:g
        [~,~,sgm_k1(:,k)] = Disc2_value(D(:,groups{k}),'',flag);
    end
elseif strcmp(flag,'WD')
    [Disc,sigma] = Disc2_value(D,q,flag);
end
sgm_k = zeros(N,N,g);
for k = 1:g
    [~,sgm_k(:,:,k)] = Disc2_value(D(:,groups{k}),'',flag);
end


D0 = D;    
Disc0 = Disc;
Disc0vec = zeros(Reps,1);
if g==1
    return;
end


% ���ȼ�����ֵ�仯�ı���,��ÿ InIter �ε�����
% ����һ����ֵ T = T*ratio
if T0 ~= 0
    ratio = (T1/T0)^(1/(ceil(AllIter/InIter)-1));
else
    ratio = 0;
end
%�����㷨
for rep = 1:Reps
    %ÿ�ζ�����һ���ظ��ѵõ����Ž���Ϊ��ʼ��,����TA�㷨
    Disc0vec(rep) = Disc;
    T = T0;
    for Oiter = 1:AllIter
        if mod(Oiter,InIter)==0
            T = T*ratio;
        end
        %����Ҫ�����е�group���k, �Լ�Ҫ����������i1,i2
        k = randi(g-1)+1; %��һ��group���ý������û�
        i = randi(N);
        j = randi(N);
        while j == i
            j=randi(N);
        end
        
        %�������ڼ��� delta ��alpha1,alpha2,beta
        if strcmp(flag,'CD') || strcmp(flag,'MD')
            alpha1 = sgm_k1(j,k)/sgm_k1(i,k);
            alpha2 = sgm_k(j,j,k)/sgm_k(i,i,k);
        end
        beta = sgm_k(:,j,k)./sgm_k(:,i,k);
        
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
            delta = delta+(alpha1-1)*(sigma1(j)/alpha1-sigma1(i))...
                +(alpha2-1)*(sigma(i,i)-sigma(j,j)/alpha2);
        end
        
        %�ж��Ƿ����
        u = rand();
        if delta <  Disc*T-epsilon %|| u < 1e-3
            %D_ = D;
            %[Disc_,sigma_,sigma1_] = Disc2_value(D_,'',flag);
            % ���� D,Disc
            temp=D(i,groups{k});
            D(i,groups{k})=D(j,groups{k});
            D(j,groups{k})=temp;
            Disc = Disc + delta;
            % ���� sigma1,diag(sigma),sgm_k1,diag(sgm_k)
            if strcmp(flag,'CD') || strcmp(flag,'MD')
                sigma1(i)=sigma1(i)*alpha1;
                sigma1(j)=sigma1(j)/alpha1;
                sigma(i,i)=sigma(i,i)*alpha2;
                sigma(j,j)=sigma(j,j)/alpha2;
                temp = sgm_k1(i,k);
                sgm_k1(i,k) = sgm_k1(j,k);
                sgm_k1(j,k) = temp;
                temp = sgm_k(i,i,k);
                sgm_k(i,i,k) = sgm_k(j,j,k);
                sgm_k(j,j,k) = temp;
            end
            % ���� sigma �� sgm_k �ǶԽǲ���
            for t=1:N
                if t~=i && t~=j
                    sigma(i,t)=sigma(i,t)*beta(t);
                    sigma(j,t)=sigma(j,t)/beta(t);
                    sigma(t,i)=sigma(i,t);
                    sigma(t,j)=sigma(j,t);
                    temp = sgm_k(i,t,k);
                    sgm_k(i,t,k) = sgm_k(j,t,k);
                    sgm_k(j,t,k) = temp;
                    sgm_k(t,i,k) = sgm_k(i,t,k);
                    sgm_k(t,j,k) = sgm_k(j,t,k);
                end
            end
            
            %delta_ = WD2_value(D)-WD2_value(D_);
            %if abs(delta-delta_)>epsilon
            %    fprintf('Wrong delta!\n');
            %end
            
            
            % �� delta<0 ˵����ʼ�½�.
            if delta < -epsilon
                % �������ȱ��θ��Ž⣬���¼֮
                if Disc < Disc0vec(rep) - epsilon
                    Disc0vec(rep) = Disc;
                    % ������ȫ�ָ��Ž⣬���¼֮
                    if Disc < Disc0-epsilon
                        Disc0 = Disc;
                        D0 = D;
                    end
                end
            end
        end
        if iswrite && rep == rep0
            fprintf(fid,'%.8f\n',Disc);
        end
    end
end

if iswrite
    fclose(fid);
end

end



%{
%���Դ���
%N = 16; q = [4;4;4;4;4];
%N = 36; q = [3;3;3;3;6];
N = 25; q = [5;5;5;5;5;5];
%N = 32; q = [4;4;4;4;4;4];
%N = 50; q = 5*ones(6,1);
InIter = 1000;
OutIter = 100*InIter;
T0 = 0.01; T1 = 1e-6;
Reps = 10;
[Disc0, Disc_vec, D0] = TA_MAD(N,q,OutIter,InIter,T0,T1,30,'out.txt');
Disc = Disc_LevelPerm(D0,q);
fprintf('%.6f  %.6f  %.6f  %d\n', Disc0, Disc, mean(Disc_vec), isOA(D0,2));
Disc_seq = importdata('out.txt');
plot(1:length(Disc_seq),Disc_seq);
%}



