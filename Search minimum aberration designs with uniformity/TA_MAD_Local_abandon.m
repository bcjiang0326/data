function [aveCD1, aveCD_vec, D1] = TA_MAD_Local_abandon(D0,Jset,OutIter,InIter,T0,T1,Reps,filename)
%20181225 by Bochuan Jiang
%20190615 ����
%�ο� Jiang and Ai (2019) Construction of uniform asymmetric designs
%����TA�ֲ������㷨
%      ��ʼ��� D0 �����谴ˮƽ�������У�����ͬˮƽ�������������С�
%      �ڵ��������н����Ըı� Jset ��Ԫ�ر�ǵ��У���˳�Ϊ�ֲ�������
%ע����Ӧ�ýǶ�������
%   (1) TA_MAD_Local ����������ȫ��� TA_MAD_Recur ���������������TA_MAD_Local 
%       �ⲿ�����D0�����У��������� TA_MAD_Recur �Ĵ������о��������б�ָ��Ϊ�ɱ���Jset��
%       �������ĺô��ǿ��Կ����������������������һ�ν�������һ�е� TA_MAD_Recur��
%   (2) ��Jset=[1,...,n]'ʱ��TA_MAD_Local ��ͬ�� TA_MAD��
%INPUT:
%   D0: N-by-n matrix һ����ʼ���
%   Jset: m-by-1 vector, m <= n, Ԫ���������У����D0�пɸı����
%   OutIter: �ⲿ��������
%   InIter: �ڲ���������
%   T0: ��ʼ��ֵ��λ��(0,1)֮�䣬����1e-2
%   T1: ������ֵ��λ��(0,1)֮�䣬����1e-5
%   Reps: �ظ��������� (default,1)
%   filename: �ṩ�ļ����Ļ����򽫽��д���ļ�
%OUTPUT:
%   aveCD1: minimum average discrepancy
%   aveCD_vec: Reps�ظ������õ��� minimum average discrepancy ����
%   D1: MAD design

[N,n] = size(D0);
q = zeros(n,1);
for i = 1:n
    q(i) = length(unique(D0(:,i)));
end
[isgrouped,uq,m,group_id] = isgrouped_ByLevNum(q);
if ~isgrouped
    error('TA_MAD_Local: q is not grouped by number of levels!\n');
end
if ~issorted(Jset) || any(Jset~=floor(Jset)) || max(Jset) > n || min(Jset) < 1 ...
        || length(unique(Jset))~=length(Jset) || size(Jset,2) ~= 1
    error('TA_MAD_Local: Wrong Jset!\n');
end
if OutIter < 1 || InIter < 1
    error('TA_MAD_Local: Wrong OutIter, InIter!\n');
end
if ~(0 <= T1 && T1 <= T0 && T0 < 1)
    error('TA_MAD_Local: Wrong T0, T1!\n');
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

aveCD1 = inf;
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
D = D0;
for rep = 1:Reps
    % ���ȶ�Jset��ǵ��н������������Ϊ��rep�ظ��ĳ�ʼ����
    for j = Jset'
        D(:,j) = D0(randperm(N,N),j);
    end
    % ���ɳ�ʼ��sigma
    sigma = zeros(N,N);
    for i = 1:N-1
        for i2 = i+1:N
            vec = abs(D(i,:)-D(i2,:)) < epsilon;
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
    if aveCD_vec(rep) < aveCD1-epsilon
        aveCD1 = aveCD_vec(rep);
        if nargout > 2
            D1=D;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %����TA�㷨
    T = T0;
    for Oiter = 1:OutIter
        if mod(Oiter,InIter)==0
            T = T*ratio;
        end
        %����Ҫ��������j����������Ԫ��(i,j),(i2,j)
        j = Jset(randi(length(Jset)));
        i = randi(N);
        i2 = randi(N);
        while i2 == i || D(i,j)==D(i2,j)
            i2=randi(N);
        end
        %�����������ˮƽ��
        k = group_id(j);
        
        %����delta
        delta = 0;
        for t = 1:N
            if t ~= i && D(t,j)==D(i,j)
                delta = delta + sigma(i2,t)-sigma(i,t)/c2(k);
            elseif t~=i2 && D(t,j)==D(i2,j)
                delta = delta + sigma(i,t)-sigma(i2,t)/c2(k);
            end
        end
        delta = delta*2*(c2(k)-1);
        
        %�ж��Ƿ����
        if delta <  aveCD*T-epsilon
            % ���� sigma
            for t = 1:N
                if t ~= i && D(t,j) == D(i,j)
                    sigma(i,t) = sigma(i,t)/c2(k);
                    sigma(i2,t) = sigma(i2,t)*c2(k);
                    sigma(t,i)=sigma(i,t);
                    sigma(t,i2)=sigma(i2,t);
                elseif t ~= i2 && D(t,j) == D(i2,j)
                    sigma(i,t) = sigma(i,t)*c2(k);
                    sigma(i2,t) = sigma(i2,t)/c2(k);
                    sigma(t,i)=sigma(i,t);
                    sigma(t,i2)=sigma(i2,t);
                end
            end
            % ���� D,aveCD
            temp=D(i,j);D(i,j)=D(i2,j);D(i2,j)=temp;
            aveCD = aveCD + delta;
            % �� delta<0 ˵����ʼ�½�.
            if delta < -epsilon
                % �������ȱ��θ��Ž⣬���¼֮
                if aveCD < aveCD_vec(rep) - epsilon
                    aveCD_vec(rep) = aveCD;
                    % ������ȫ�ָ��Ž⣬���¼֮
                    if aveCD_vec(rep) < aveCD1-epsilon
                        aveCD1 = aveCD_vec(rep);
                        if nargout > 2
                            D1 = D;
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

k = 2; m = 2; p = 2; 
s = p^m; 
prim_poly = gfprimfd(m,'min',p);
[D0,GenMat,OA_polys,GM_polys] = RH_OA_pow(k,m,p,prim_poly);
[N,n] = size(D0);
q = s*ones(n,1);
Jset = [n-1;n];

InIter = 1000;
OutIter = 100*InIter;
T0 = 0.01; T1 = 1e-6;
Reps = 10;
filename = 'out.txt';

[aveCD1, aveCD_vec, D1] = TA_MAD_Local(D0,Jset,OutIter,InIter,T0,T1,Reps,'out.txt');
aveCD = aveCD_LevelPerm(D1,q);
fprintf('%.6f  %.6f  %.6f  %d\n', aveCD1, aveCD, mean(aveCD_vec), isOA(D1,2));
aveCD_seq = importdata('out.txt');
plot(1:length(aveCD_seq),aveCD_seq);
%}



