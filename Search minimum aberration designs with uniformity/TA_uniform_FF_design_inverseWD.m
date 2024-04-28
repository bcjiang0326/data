function [WD0, WDvec, D0] = TA_uniform_FF_design_inverseWD(D,AllIter,InIter,T0,T1,Reps,writedata)
% 20200213 ��WD������ƣ�����ο����� TA_uniform_FF_design(...)      .
% INPUT:
%       D: N-by-n balanced design
%       AllIter: �ܵ�������=OutIter*InIter
%       InIter: �ڲ���������
%       T0: ��ʼ��ֵ��λ��(0,1)֮�䣬����1e-2
%       T1: ������ֵ��λ��(0,1)֮�䣬����1e-6
%       Reps: �ظ����� (default,1)
%       writedata: �Ƿ񽫹�������д���ļ��� 1,д��; 0(default) ��д
% OUTPUT:
%       WD0: �õ����ž���� Squared WD
%       WDvec: Reps-by-1 vector, ��¼ÿ���ظ������Ž�
%       D0: �õ������ž���

[N,n] = size(D);
[isbal,q] = isBalanced(D);
if ~isbal
    error('D is not balanced!\n');
end
if nargin < 6
    Reps = 1;
    writedata = 0;
elseif nargin < 7
    writedata = 0;
end
if Reps > 1
    writedata = 0;
end

% �����ļ����
if writedata
    fid = fopen('T_Disc_Disc.txt','w');
end

% ���ù������
epsilon = 1e-12;

% ���ȼ�����ֵ�仯�ı���,��ÿ InIter �ε�����
% ����һ����ֵ T = T*ratio
ratio = (T1/T0)^(1/(ceil(AllIter/InIter)-1));

KK = N./q; % ��¼�����У�ÿ��ˮƽ�ظ�����

D0 = D;
WD0 = WD2_value(D0,q);
WDvec = zeros(Reps,1);

for rep = 1:Reps
    % level_id{j}��¼��j���и�Ԫ�س��ֵ�λ��
    level_id = cell(1,n);
    for j = 1:n
        [~,id] = sort(D(:,j));
        level_id{j} = reshape(id,[N/q(j),q(j)]);
    end
    % �������ɳ�ʼ�� sigma
    sigma = zeros(N,N);
    for i = 1:N-1
        for j = i+1:N
            avec = abs(D(i,:)-D(j,:))./q';
            sigma(i,j) = 1/N^2*prod(1.5-avec.*(1-avec));
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
    WDvec(rep) = Disc;
    if WDvec(rep) > WD0+epsilon
        WD0 = WDvec(rep);
        if nargout > 2
            D0=D;
        end
    end
    
    %�����㷨
    T = T0;
    for Oiter = 1:AllIter
        if mod(Oiter,InIter)==0
            T = T*ratio;
        end
        %����Ҫ��������
        k=randi(n);
        %k = n;
        %����Ҫ����������ˮƽlv1,lv2
        lv1=randi(q(k))-1;
        lv2=randi(q(k))-1;
        while lv2==lv1
            lv2=randi(q(k))-1;
        end
        if lv1 > lv2
            temp = lv1; lv1 = lv2; lv2 = temp;
        end
        %�û�2ˮƽ���ӡ�3ˮƽ���ӡ�����4ˮƽ���ӵ�0��2������4ˮƽ���ӵ�1��3�������ı�WDֵ.
        if q(k)==2 || q(k)==3 || ( q(k) == 4 && lv1 == 0 && lv2 == 2 ) || ( q(k) == 4 && lv1 == 1 && lv2 == 3 )
            if writedata
                fprintf(fid,'%.8f %.8f %.8f %.8f %.8f\n',0,T,Disc,WD0,Oiter);
            end
            continue;
        end
        Iset = level_id{k}(:,lv1+1);
        Jset = level_id{k}(:,lv2+1);
        Tset = reshape(level_id{k}(:,[1:lv1,lv1+2:lv2,lv2+2:end]),[N-2*KK(k),1]);
        
        %�������ڼ��� delta �� beta
        ivec = abs(lv1-D(Tset,k))/q(k);
        jvec = abs(lv2-D(Tset,k))/q(k);
        beta = (1.5 - jvec.*(1-jvec))./(1.5-ivec.*(1-ivec));
        
        % ���� delta
        delta = ( sum(sigma(Iset,Tset)) - sum(sigma(Jset,Tset))./beta')*(beta-1)*2;
        
        %�ж��Ƿ����
        uu = rand();
        if delta >  -Disc*T*(1+epsilon) || uu < 1e-3
            % ���� D,Disc,level_id
            D(Iset,k) = lv2; D(Jset,k) = lv1;
            Disc = Disc + delta;
            level_id{k}(:,lv1+1) = Jset;
            level_id{k}(:,lv2+1) = Iset;
            % ���� sigma
            sigma(Tset,Iset) = sigma(Tset,Iset).*(beta*ones(1,KK(k)));
            sigma(Iset,Tset) = sigma(Tset,Iset)';
            sigma(Tset,Jset) = sigma(Tset,Jset)./(beta*ones(1,KK(k)));
            sigma(Jset,Tset) = sigma(Tset,Jset)';
            % �� delta>0 ˵����ʼ����.
            if delta > Disc*epsilon
                % �������ȱ��θ���⣬���¼֮
                if Disc > WDvec(rep)*(1 + epsilon)
                    WDvec(rep) = Disc;
                    % ������ȫ�ָ��Ž⣬���¼֮
                    if WDvec(rep) > WD0*(1 + epsilon)
                        WD0 = WDvec(rep);
                        D0 = D;
                        fprintf('inv LP: N=%d, n=%d, rep=%d, iter=%d: %.6f\n',N,n,rep,Oiter,WD0);
                    end
                end
            end
        end
        if writedata
            fprintf(fid,'%.8f %.8f %.8f %.8f %.8f\n',delta/Disc,T,Disc,WD0,Oiter);
        end
    end
    %fprintf('LP: N = %d, n = %d, rep = %d: ******** %.6f, %.6f\n',N,n,rep,WDvec(rep),WD0);
end

if writedata
    fclose(fid);
end

end

%{
oa = importdata('OA from web/MA.36.3.7.6.3.finney.txt');
oa = oa(:,5:10);
aveWD = aveDisc_LevelPerm(oa,'','WD');

InIter = 100;
OutIter = 100*InIter;

[WD0, WDvec, D0] = TA_uniform_FF_design_inverseWD(oa,OutIter,InIter,1e-3,1e-6,1,1);
D0 = sortrows(D0);
fprintf('%.8f\n',WDvec);
fprintf('%.8f\n',WD2_value(D0));
out = importdata('T_Disc_Disc.txt');
plot(out(:,5),out(:,3));
%}