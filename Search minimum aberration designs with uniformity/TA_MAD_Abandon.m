function [aveCD0, aveCD_vec, D0] = TA_MAD_Abandon(N,q,OutIter,InIter,T0,T1,Reps,filename)
%20181211 by Bochuan Jiang
%20181217 ע������ʵ��Ч�ʲ���TA_MAD, ������á�
%�ο� Jiang and Ai (2019) Construction of uniform asymmetric designs
%����TA�㷨����
%      ��� D �ĸ����Ӳ��ǰ���ˮƽ���������У����������� D �ĸ��У�
%      ʹ�ø���ˮƽ��Ϊ��������
%INPUT:
%   N: �������
%   q: n-by-1 vector, ��Ƹ�����ˮƽ��
%   OutIter: �ⲿ��������
%   InIter: �ڲ���������
%   T0: ��ʼ��ֵ��λ��(0,1)֮�䣬����1e-2
%   T1: ������ֵ��λ��(0,1)֮�䣬����1e-5
%   Reps: �ظ��������� (default,1)
%   filename: �ṩ�ļ����Ļ����򽫽��д���ļ�
%OUTPUT:
%   aveCD0: minimum average discrepancy
%   aveCD_vec: Reps�ظ������õ��� minimum average discrepancy ����
%   D0: MAD design

N = floor(N);
if N < 1
    error('TA_MAD: Wrong N!\n');
end
[isgrouped,uq,m,group_id] = isgrouped_ByLevNum(q);
if ~isgrouped
    error('TA_MAD: q is not grouped by number of levels!\n');
end
if OutIter < 1 || InIter < 1
    error('TA_MAD: Wrong OutIter, InIter!\n');
end
if ~(0 <= T1 && T1 <= T0 && T0 < 1)
    error('TA_MAD: Wrong T0, T1!\n');
end
if nargin < 7 || Reps < 1
    Reps = 1;
end
if nargin < 8
    iswrite = 0;
else
    iswrite = 1; % �������
    fid = fopen(filename,'w');
    %��¼ Reps �ظ��л����СaveCD ����������
    aveCD_seq0 = zeros(OutIter,1);
    aveCD_seq = zeros(OutIter,1);
    rep0 = 1; %��¼���Ž��Ӧ��rep
end

n = length(q);
epsilon = 1e-11;

aveCD0 = inf;
aveCD_vec = zeros(Reps,1);

alpha1 = 2; %������
for k = 1:length(uq)
    if mod(uq(k),2)==0
        alpha1 = alpha1*((26*uq(k)^2+1)/24/uq(k)^2)^m(k);
    else
        alpha1 = alpha1*((13*uq(k)^2-1)/12/uq(k)^2)^m(k);
    end
end
alpha1 = (13/12)^n-alpha1;

c1 = zeros(size(uq));
for k = 1:length(uq)
    if mod(uq(k),2)==0
        c1(k) = (13*uq(k)-2)*(uq(k)-1)/12;
    else
        c1(k) = (13*uq(k)^2-2*uq(k)-3)*(uq(k)-1)/12/uq(k);
    end
end
alpha2 = prod((c1./uq./(uq-1)).^m)/N^2; %���������ϵ��

c2 = zeros(size(uq));
for k = 1:length(uq)
    if mod(uq(k),2)==0
        c2(k) = 15*uq(k)/(13*uq(k)-2);
    else
        c2(k) = (15*uq(k)^2-3)/(13*uq(k)^2-2*uq(k)-3);
    end
end



% ���ȼ�����ֵ�仯�ı���,��ÿ InIter �ε�����
% ����һ����ֵ T = T*ratio
if T0 ~= 0
    ratio = (T1/T0)^(1/(ceil(OutIter/InIter)-1));
else
    ratio = 0;
end
%�����㷨
for rep = 1:Reps
    % ��������һ��������ظ��ľ��� D
    D = rand_U_Type_design(q,N);
    % ���� D �и�Ԫ�أ���ˮƽ��Ӧ���к�
    RowsPerLevel = cell(1,n);% RowsPerLevel{j}(l+1,:) ��ʾ��j�����ӵ�lˮƽ�������к�
    for j = 1:n
        RowsPerLevel{j} = zeros(q(j),N/q(j));
        for l = 0:q(j)-1 %ע�⣬ˮƽ�Ǵ�0��ʼ�����½Ǳ��1��ʼ
            RowsPerLevel{j}(l+1,:) = find(D(:,j)==l)';
        end
    end
    % �����ʼ���� D �� sigma
    sigma = zeros(N,N);
    for i = 1:N-1
        for i2 = i+1:N
            vec = abs(D(i,:)-D(i2,:))<epsilon;
            v = zeros(length(m),1);
            %v ��¼ i �к� i2 �и���group�� Hamming distance,����ͬλ�õĸ���
            bg = 1; ed = m(1);
            v(1) = sum(vec(bg:ed));
            for k = 2:length(m)
                bg = bg+m(k-1);
                ed = bg+m(k)-1;
                v(k) = sum(vec(bg:ed));
            end
            sigma(i,i2) = alpha2*prod(c2.^v);
            sigma(i2,i) = sigma(i,i2);
        end
    end
    for i = 1:N
        sigma(i,i) = alpha2*prod(c2.^m);
    end
    
    % ��ʼ��aveCD
    aveCD=0;
    for i = 1:N-1
        for i2 = i+1:N
            aveCD = aveCD + sigma(i,i2);
        end
    end
    aveCD = aveCD*2;
    for i = 1:N
        aveCD = aveCD + sigma(i,i);
    end
    aveCD = alpha1+aveCD;
    
    % �жϵ����ظ��ĳ�ʼ���Ƿ����
    aveCD_vec(rep) = aveCD;
    if aveCD_vec(rep) < aveCD0-epsilon
        aveCD0 = aveCD_vec(rep);
        if nargout > 2
            D0=D;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %����TA�㷨
    T = T0;
    for Oiter = 1:OutIter
        if mod(Oiter,InIter)==0
            T = T*ratio;
        end
        %����Ҫ�������� j0
        j0 = randi(n);
        %�����������ˮƽ��
        k = group_id(j0);
        %��������������Ԫ��(i1,j0),(i2,j0)
        l1 = randi(q(j0))-1;%�����������������ͬˮƽl1,l2
        l2 = randi(q(j0))-1;
        while l2 == l1
            l2 = randi(q(j0))-1;
        end
        r1 = randi(N/q(j0));%l1ˮƽ�������еĵڼ���
        r2 = randi(N/q(j0));%l2ˮƽ�������еĵڼ���
        i1 = RowsPerLevel{j0}(l1+1,r1);
        i2 = RowsPerLevel{j0}(l2+1,r2);
        
        %����delta
        delta = 0;
        for t = RowsPerLevel{j0}(l1+1,:)
            if t ~= i1 
                delta = delta + sigma(i2,t)-sigma(i1,t)/c2(k);
            end
        end
        for t = RowsPerLevel{j0}(l2+1,:)
            if t~=i2 
                delta = delta + sigma(i1,t)-sigma(i2,t)/c2(k);
            end
        end
        delta = delta*2*(c2(k)-1);
        
        %�ж��Ƿ����
        if delta <  aveCD*T-epsilon
            % ���� sigma
            for t = RowsPerLevel{j0}(l1+1,:)
                if t ~= i1
                    sigma(i1,t) = sigma(i1,t)/c2(k);
                    sigma(i2,t) = sigma(i2,t)*c2(k);
                    sigma(t,i1) = sigma(i1,t);
                    sigma(t,i2) = sigma(i2,t);
                end
            end
            for t = RowsPerLevel{j0}(l2+1,:)
                if t~=i2
                    sigma(i1,t) = sigma(i1,t)*c2(k);
                    sigma(i2,t) = sigma(i2,t)/c2(k);
                    sigma(t,i1)=sigma(i1,t);
                    sigma(t,i2)=sigma(i2,t);
                end
            end
            % ���� D,aveCD
            temp=D(i1,j0);D(i1,j0)=D(i2,j0);D(i2,j0)=temp;
            aveCD = aveCD + delta;
            % ����RowsPerLevel
            RowsPerLevel{j0}(l1+1,r1) = i2;
            RowsPerLevel{j0}(l2+1,r2) = i1;
            
            % �� delta<0 ˵����ʼ�½�.
            if delta < -epsilon
                % �������ȱ��θ��Ž⣬���¼֮
                if aveCD < aveCD_vec(rep) - epsilon
                    aveCD_vec(rep) = aveCD;
                    % ������ȫ�ָ��Ž⣬���¼֮
                    if aveCD_vec(rep) < aveCD0-epsilon
                        aveCD0 = aveCD_vec(rep);
                        if nargout > 2
                            D0 = D;
                        end
                        if iswrite
                            rep0 = rep;
                        end
                    end
                end
            end
        end
        if iswrite 
            aveCD_seq(Oiter) = aveCD;
        end
    end
    if iswrite && rep0 == rep
        aveCD_seq0 = aveCD_seq;
    end
end

if iswrite
    fprintf(fid,'%.8f\n',aveCD_seq0);
    fclose(fid);
end

end



%{
%���Դ���
%N = 16; q = [4;4;4;4;4];
N = 36; q = [3;3;3;3;6];
%N = 25; q = [5;5;5;5;5;5];
%N = 32; q = [4;4;4;4;4;4];
%N = 50; q = 5*ones(6,1);
InIter = 1000;
OutIter = 100*InIter;
T0 = 0.01; T1 = 1e-6;
%T0 = 0; T1 = 0;
Reps = 10;
[aveCD0, aveCD_vec, D0] = TA_MAD(N,q,OutIter,InIter,T0,T1,30,'out.txt');
aveCD = aveCD_LevelPerm(D0,q);
fprintf('%.6f  %.6f  %.6f  %d\n', aveCD0, aveCD, mean(aveCD_vec), isOA(D0,2));
aveCD_seq = importdata('out.txt');
plot(1:length(aveCD_seq),aveCD_seq);
%}



