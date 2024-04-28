function [aveDisc2, aveDisc2_vec, PI2] = TA_MAD_Combin_abandon(D0,D1,OutIter,InIter,T0,T1,flag,Reps,filename)
%20190615 by Bochuan Jiang
%�ο� Jiang and Ai (2019) Construction of uniform minimum aberration designs
%����TA�ֲ������㷨
%      D0�̶����� D2 �������û���ʹ�� D2= [D0,D2] ��average discrepancy ��С��
%INPUT:
%   D0: N-by-n0 matrix 
%   D1: N-by-n1 matrix
%   OutIter: �ⲿ��������
%   InIter: �ڲ���������
%   T0: ��ʼ��ֵ��λ��(0,1)֮�䣬����1e-2
%   T1: ������ֵ��λ��(0,1)֮�䣬����1e-5
%   flag: 'CD'(default), 'WD' or 'MD'
%   Reps: �ظ��������� (default,1)
%   filename: �ṩ�ļ����Ļ����򽫽��д���ļ�
%OUTPUT:
%   aveDisc2: minimum average discrepancy
%   aveDisc2_vec: Reps�ظ������õ��� minimum average discrepancy ����
%   PI2: ���ŵ��û���D2 = [D0,D1(PI2,:)]
epsilon = 1e-10;
[N,n0] = size(D0);
[N1,n1] = size(D1);
if N ~= N1
    error('D0 and D1 have distinct numbers of rows!\n');
end
q0 = zeros(n0,1);
for i = 1:n0
    q0(i) = length(unique(D0(:,i)));
end
q1 = zeros(n1,1);
for i = 1:n1
    q1(i) = length(unique(D1(:,i)));
end
if OutIter < 1 || InIter < 1
    error('TA_MAD_Combin: Wrong OutIter, InIter!\n');
end
if ~(0 <= T1 && T1 <= T0 && T0 < 1)
    error('TA_MAD_Combin: Wrong T0, T1!\n');
end
if nargin < 7
    flag = 'CD';
end
if ~strcmp(flag,'CD') && ~strcmp(flag,'WD') && ~strcmp(flag,'MD')
    error('Wrong flag!\n');
end
if nargin < 8 || Reps < 1
    Reps = 1;
end
if nargin < 9
    iswrite = 0;
else
    iswrite = 1; % �������
    fid = fopen(filename,'w');
    %��¼ Reps �ظ��л����СaveCD ����������
    aveDisc_seq0 = zeros(OutIter,1);
    aveDisc_seq = zeros(OutIter,1);
    rep0 = 1; %��¼���Ž��Ӧ��rep
end

% ���ȼ�����ֵ�仯�ı���,��ÿ InIter �ε�����
% ����һ����ֵ T = T*ratio
if T0 ~= 0
    ratio = (T1/T0)^(1/(ceil(OutIter/InIter)-1));
else
    ratio = 0;
end
%�����㷨
q = [q0;q1];
%ȫ�����Ž� aveDisc2,PI2
aveDisc2 = aveDisc_LevelPerm([D0,D1],q,flag);
aveDisc2_vec = zeros(Reps,1);
PI2 = (1:N)';
for rep = 1:Reps    
    %������rep�ε����ĳ�ʼ��(starting design)��
    %��Ϊѭ������Ҫ�ı����D1_loop, PI_loop, aveDisc_loop
    D1_loop = D1;
    PI_loop = (1:N)';
    aveDisc_loop = aveDisc_LevelPerm([D0,D1_loop],q,flag); 
    %��rep�����У��ֲ����ŵļ�¼����aveDisc_vec(rep),PI_opt
    aveDisc2_vec(rep) = aveDisc_loop;
    PI_opt = PI_loop;
    %����ֲ����Ž�����ȫ�����Ž⣬����ȫ�����Ž�
    if aveDisc2_vec(rep) < aveDisc2
        aveDisc2 = aveDisc2_vec(rep);
        PI2 = PI_opt;
    end
    
    %����TA�㷨
    T = T0;
    for Oiter = 1:OutIter
        if mod(Oiter,InIter)==0
            T = T*ratio;
        end
        %����D1_loop��Ҫ������ i,j ��
        i = randi(N);
        j = randi(N);
        while j == i
            j=randi(N);
        end
        temp = D1_loop(i,:);
        D1_loop(i,:) = D1_loop(j,:);
        D1_loop(j,:) = temp;
        
        %���㽻�����к��µ� aveDisc_tilde
        aveDisc_tilde = aveDisc_LevelPerm([D0,D1_loop],q,flag);
        delta = aveDisc_tilde-aveDisc_loop;
        
        %�ж��Ƿ����
        if delta >= aveDisc_loop*T-epsilon %�����£��˻�
            temp = D1_loop(i,:);
            D1_loop(i,:) = D1_loop(j,:);
            D1_loop(j,:) = temp;
        else % ���� PI_loop,aveDisc_loop
            temp = PI_loop(i); PI_loop(i) = PI_loop(j); PI_loop(j) = temp;
            aveDisc_loop = aveDisc_tilde;
            
            % �� delta<0 ˵����ʼ�½�.
            if delta < -epsilon
                % �������ȱ��θ��Ž⣬���¼֮
                if aveDisc_loop < aveDisc2_vec(rep) - epsilon
                    aveDisc2_vec(rep) = aveDisc_loop;
                    PI_opt = PI_loop;
                    % ������ȫ�ָ��Ž⣬���¼֮
                    if aveDisc2_vec(rep) < aveDisc2-epsilon
                        aveDisc2 = aveDisc2_vec(rep);
                        PI2 = PI_opt;
                        fprintf('%d  %d  %.6f\n',rep,Oiter,aveDisc2);                       
                        if iswrite
                            rep0 = rep;
                        end
                    end
                end
            end
        end
        if iswrite 
            aveDisc_seq(Oiter) = aveDisc_loop;
        end
    end
    if iswrite && rep0 == rep
        aveDisc_seq0 = aveDisc_seq;
    end
end

if iswrite
    fprintf(fid,'%.8f\n',aveDisc_seq0);
    fclose(fid);
end

end



%{
%���Դ���
InIter = 1000;
OutIter = 100*InIter;
T0 = 0.01; T1 = 1e-6;
flag = 'WD';
Reps = 1;
filename = 'out.txt';

D0 = importdata('OA from web/oa.32.9.4.2.a.txt');
D0 = importdata('OA from web/oa.18.7.3.2.txt');
D1 = D0;

[aveDisc2, aveDisc2_vec, PI2] = TA_MAD_Combin(D0,D1,OutIter,InIter,T0,T1,flag,Reps,'out.txt');
D2 = [D0,D1(PI2,:)];
aveDisc = aveDisc_LevelPerm(D2,'',flag);
fprintf('%.6f  %.6f  %.6f\n', aveDisc2, aveDisc, mean(aveDisc2_vec));
aveDisc_seq = importdata('out.txt');
plot(1:length(aveDisc_seq),aveDisc_seq);
%}



