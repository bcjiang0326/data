function [aveCD1, aveCD_vec, D1] = TA_MAD_Recur(D0,s,OutIter,InIter,T0,T1,Reps,filename)
%20181216 by Bochuan Jiang
%D0 Ϊһ��Ļ�ˮƽ��ƣ�Ҫ����ͬˮƽ������������һ�𡣱�������Ŀ���ǽ�D0����һ��sˮƽ����
%   ����TA�㷨���Ż��������У�ʹ������ƽ����������С��
%   ��D1 = [D0,dj], D0���ֲ������Ż�dj��djΪsˮƽ
%ע��TA_MAD_Recur �ɱ� TA_MAD_Local ��ȫȡ��������ʹ��TA_MAD_Local
%INPUT:
%   D0: һ���ˮƽ���
%   s: �����ӵ��е�ˮƽ��
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

[N,n0] = size(D0);
if mod(N,s)~=0
    error('TA_MAD_Recur: a balanced s-level factor can not be added to D0!\n');
end
q0 = zeros(n0,1);
for j = 1:n0
    q0(j) = length(unique(D0(:,j)));
end
[isgrouped,uq0,m0] = isgrouped_ByLevNum(q0);
if ~isgrouped
    error('TA_MAD_Recur: D0 is not grouped by number of levels!\n');
end
if OutIter < 1 || InIter < 1
    error('TA_MAD_Recur: Wrong OutIter, InIter!\n');
end
if ~(0 <= T1 && T1 <= T0 && T0 < 1)
    error('TA_MAD_Recur: Wrong T0, T1!\n');
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
    aveCD_seq1 = zeros(OutIter,1);
    aveCD_seq = zeros(OutIter,1);
    rep0 = 1; %��¼���Ž��Ӧ��rep
end

n = n0+1; %������Ƶ�����
K = find(uq0==s,1);%����ӵ�sˮƽ�����ڵ�ˮƽ��
if isempty(K)
    %�������е�ˮƽ������ӵ����һ��
    uq = [uq0;s];
    m = [m0;1];
    J = n;
    K = length(m);
else %���е�ˮƽ��, ����и�ˮƽ����������
    J = sum(m0(1:K))+1;
    m = m0; m(K) = m0(K)+1;
    uq = uq0;
end
dJ = rand_U_Type_design(s,N);
D = [D0(:,1:J-1),dJ,D0(:,J:end)];
D1 = D;
%q = cat(1,q0,q0(end));
%uq = unique(q); %��Ƴ��ֵĲ�ͬˮƽ��ˮƽ������
%m = histc(q,[uq-0.5;inf]); %��ˮƽ���Ӹ�������
%m(end) = []; %��ˮƽ���ӵĸ�����m(i) Ϊˮƽ��Ϊ uq(i) �����Ӹ���


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
for rep = 1:Reps
    % ��������һ������������ɳ�ʼ�� sigma
    D(:,J) = D(randperm(N),J);
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
    if aveCD_vec(rep) < aveCD1-epsilon
        aveCD1 = aveCD_vec(rep);
        if nargout > 2
            D1(:,J)=D(:,J);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %����TA�㷨
    T = T0;
    for Oiter = 1:OutIter
        if mod(Oiter,InIter)==0
            T = T*ratio;
        end
        %����Ҫ��������J�е�����Ԫ��(i,J),(i2,J)
        i = randi(N);
        i2 = randi(N);
        while i2 == i || D(i,J)==D(i2,J)
            i2=randi(N);
        end
        
        %����delta
        delta = 0;
        for t = 1:N
            if t ~= i && D(t,J)==D(i,J)
                delta = delta + sigma(i2,t)-sigma(i,t)/c2(K);
            elseif t ~= i2 && D(t,J)==D(i2,J)
                delta = delta + sigma(i,t)-sigma(i2,t)/c2(K);
            end
        end
        delta = delta*2*(c2(K)-1);
        
        %�ж��Ƿ����
        if delta <  aveCD*T-epsilon
            % ���� sigma
            for t = 1:N
                if t ~= i && D(t,J) == D(i,J)
                    sigma(i,t) = sigma(i,t)/c2(K);
                    sigma(i2,t) = sigma(i2,t)*c2(K);
                    sigma(t,i) = sigma(i,t);
                    sigma(t,i2) = sigma(i2,t);
                elseif t ~= i2 && D(t,J) == D(i2,J)
                    sigma(i,t) = sigma(i,t)*c2(K);
                    sigma(i2,t) = sigma(i2,t)/c2(K);
                    sigma(t,i) = sigma(i,t);
                    sigma(t,i2) = sigma(i2,t);
                end
            end
            % ���� D,aveCD
            temp=D(i,J);D(i,J)=D(i2,J);D(i2,J)=temp;
            aveCD = aveCD + delta;
            
            %aveCD_ = aveCD_LevelPerm(D,q);
            %fprintf('%e\n',aveCD_-aveCD);
            
            % �� delta<0 ˵����ʼ�½�.
            if delta < -epsilon
                % �������ȱ��θ��Ž⣬���¼֮
                if aveCD < aveCD_vec(rep) - epsilon
                    aveCD_vec(rep) = aveCD;
                    % ������ȫ�ָ��Ž⣬���¼֮
                    if aveCD_vec(rep) < aveCD1-epsilon
                        aveCD1 = aveCD_vec(rep);
                        if nargout > 2
                            D1(:,J) = D(:,J);
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
        aveCD_seq1 = aveCD_seq;
    end
end

if iswrite
    fprintf(fid,'%.8f\n',aveCD_seq1);
    fclose(fid);
end

end



%{
%���Դ���
N = 16; q = [4;4;4;4;4];
InIter = 1000;
OutIter = 100*InIter;
T0 = 0.01; T1 = 1e-5;
Reps = 30;
filename = 'out.txt';
[aveCD0, aveCD_vec0, D0] = TA_MAD(N,q,OutIter,InIter,T0,T1,30,filename);
aveCD = aveCD_LevelPerm(D0,q);
fprintf('%.6f  %.6f  %.6f  %d\n', aveCD0, aveCD, mean(aveCD_vec0), isOA(D0,2));
aveCD_seq = importdata('out.txt');
plot(1:length(aveCD_seq),aveCD_seq);

s = q(end);
[aveCD1, aveCD_vec1, D1] = TA_MAD_Recur(D0,s,OutIter,InIter,T0,T1,Reps,filename)
aveCD = aveCD_LevelPerm(D1,[q;s]);
fprintf('%.6f  %.6f  %.6f  %d\n', aveCD1, aveCD, mean(aveCD_vec1), isOA(D1,2));
aveCD_seq = importdata(filename);
plot(1:length(aveCD_seq),aveCD_seq);
%}



