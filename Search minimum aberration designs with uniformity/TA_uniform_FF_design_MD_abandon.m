function [MD0, MDvec, D0] = TA_uniform_FF_design_MD_abandon(D,q,OutIter,InIter,T0,T1,Reps,writedata)
% 20181219
% ���� MD ׼��
% ���� TA �㷨��� uniform fractional factorial design (UFFD)
% ��������ƣ��������� columnwise-pairwise ����, ÿ�ε�����ѡ��һ�У�ѡ������ˮƽ��i,j
%       ��ˮƽΪ i ��Ԫ�ر�Ϊ j, ��ˮƽΪ j �ı�Ϊ i.
% ������Tang, Xu and Lin (2012) AOS �����uniform fractional factorial design (UFFD). ���ڴ�
%       Tang and Xu (2012) JSPI �����һ�ֹ��� uniform FFD �������㷨�� 
% ʹ�÷�Χ�����ݸ�����D����������ˮƽ�û��µľ������
% Neighborhood: �� D ����Ϊ��D������һ������ˮƽ�û��õ�����Ƶļ���
%       .
% INPUT:
%       D: N-by-n balanced design
%       q: s-by-1 vector, ��Ƹ����ӵ�ˮƽ��
%       OutIter: �ⲿ��������
%       InIter: �ڲ���������
%       T0: ��ʼ��ֵ��λ��(0,1)֮�䣬����1e-2
%       T1: ������ֵ��λ��(0,1)֮�䣬����1e-6
%       Reps: �ظ����� (default,1)
%       writedata: �Ƿ񽫹�������д���ļ��� 1,д��; 0(default) ��д
% OUTPUT:
%       MD0: �õ����ž���� Squared Discrepancy
%       MD_vec: Reps-by-1 vector, ��¼ÿ���ظ������Ž�
%       D0: �õ������ž���

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
    fid = fopen('T_Disc_Disc.txt','w');
end

% ���ù������
epsilon = 1e-10;

% ���ȼ�����ֵ�仯�ı���,��ÿ InIter �ε�����
% ����һ����ֵ T = T*ratio     
ratio = (T1/T0)^(1/(ceil(OutIter/InIter)-1));

KK = N./q; % ��¼�����У�ÿ��ˮƽ�ظ�����

MD0 = inf; 
MDvec = zeros(Reps,1);

for rep = 1:Reps
    % ��������һ��������ظ��ľ������ɳ�ʼ�� z,sigma1,sigma2,sigma
    sigma1 = zeros(N,1);
    sigma2 = zeros(N,1);
    sigma = zeros(N,N);
    z = 0.25*( (D+0.5)./(ones(N,1)*q')-0.5);
    for i = 1:N-1
        for j = i+1:N
            sigma(i,j) = 1/N^2*prod(1.875-abs(z(i,:))-abs(z(j,:))-3*abs(z(i,:)-z(j,:))+8*(z(i,:)-z(j,:)).^2);
            sigma(j,i)=sigma(i,j);
        end 
    end
    for i = 1:N
        sigma1(i) = 2/N*prod(5/3-abs(z(i,:))-4*z(i,:).^2);
        sigma2(i) = 1/N^2*prod(1.875-2*abs(z(i,:)));
    end
    
    % ��ʼ��Disc
    Disc=0;
    for i = 1:N-1
        Disc = Disc + sum(sigma(i,i+1:N));
    end
    Disc=(19/12)^n+2*Disc+sum(sigma2)-sum(sigma1);

    % �жϵ����ظ��ĳ�ʼ���Ƿ����
    MDvec(rep) = Disc; 
    if MDvec(rep) < MD0-epsilon
        MD0 = MDvec(rep);
        D0=D;
    end
    
    %�����㷨
    T = T0;
    for Oiter = 1:OutIter
        if mod(Oiter,InIter)==0
            T = T*ratio;
        end
        %����Ҫ��������
        k=randi(n);
        %����Ҫ������ˮƽ��
        i=randi(N);
        j=randi(N);
        while j==i || D(i,k)==D(j,k)    
            j=randi(N);
        end
        %�û�2ˮƽ���ӵ�0��1����3ˮƽ���ӵ�0��2�������ı�CDֵ.
        if q(k)==2 || ( q(k)==3 && D(i,k)+D(j,k) == 2 )
            if writedata
                fprintf(fid,'%.8f %.8f %.8f %.8f %.8f\n',0,T,Disc,MD0,Oiter);
            end
            continue;
        end
        Iset = D(:,k)==D(i,k); 
        Jset = D(:,k)==D(j,k);
        Tset = ~( Iset | Jset );
      
        %�������ڼ��� delta �� alpha1,alpha2,beta
        alpha1=(5/3-abs(z(j,k))+4*z(j,k)^2)/(5/3-abs(z(i,k))+4*z(i,k)^2);
        alpha2=(1.875-2*abs(z(j,k)))/(1.875-2*abs(z(i,k)));
        beta = (1.875-abs(z(j,k))-abs(z(Tset,k))-3*abs(z(j,k)-z(Tset,k))+8*(z(j,k)-z(Tset,k)).^2)...
            ./(1.875-abs(z(i,k))-abs(z(Tset,k))-3*abs(z(i,k)-z(Tset,k))+8*(z(i,k)-z(Tset,k)).^2);
        
        % ���� delta
        delta = ( sum(sigma(Iset,Tset)) - sum(sigma(Jset,Tset))./beta')*(beta-1)*2;
        delta = delta+(alpha1-1)*(sum(sigma1(Jset))/alpha1-sum(sigma1(Iset)))...
            +(alpha2-1)*( (sum(sigma2(Iset))-sum(sigma2(Jset))/alpha2) );
        cc1 = 0; cc2 = 0;
        for lp = 1:KK(k)-1
            cc1 = cc1 + sum(diag(sigma(Iset,Iset),lp));
            cc2 = cc2 + sum(diag(sigma(Jset,Jset),lp));
        end
        delta = delta+2*(alpha2-1)*(cc1-cc2/alpha2);
        
        %�ж��Ƿ����
        if delta <  Disc*T-epsilon
            % ���� D,Disc
            temp = D(i,k); D(Iset,k) = D(j,k); D(Jset,k) = temp;
            Disc = Disc + delta;
            % ���� z,sigma1,sigma2,sigma
            temp=z(i,k);z(Iset,k)=z(j,k);z(Jset,k)=temp;            
            sigma1(Iset)=sigma1(Iset)*alpha1;
            sigma1(Jset)=sigma1(Jset)/alpha1;
            sigma2(Iset)=sigma2(Iset)*alpha2;
            sigma2(Jset)=sigma2(Jset)/alpha2;
            sigma(Tset,Iset) = sigma(Tset,Iset).*(beta*ones(1,KK(k)));
            sigma(Iset,Tset) = sigma(Tset,Iset)';
            sigma(Tset,Jset) = sigma(Tset,Jset)./(beta*ones(1,KK(k)));
            sigma(Jset,Tset) = sigma(Tset,Jset)';
            sigma(Iset,Iset) = sigma(Iset,Iset)*alpha2;
            sigma(Jset,Jset) = sigma(Jset,Jset)/alpha2;
            % �� delta<0 ˵����ʼ�½�.
            if delta < -epsilon
                % �������ȱ��θ��Ž⣬���¼֮
                if Disc < MDvec(rep) - epsilon
                    MDvec(rep) = Disc;
                    % ������ȫ�ָ��Ž⣬���¼֮
                    if MDvec(rep) < MD0-epsilon
                        MD0 = MDvec(rep);
                        if nargout > 2
                            D0 = D;
                        end
                    end
                end
            end
        end
        if writedata
            fprintf(fid,'%.8f %.8f %.8f %.8f %.8f\n',delta/Disc,T,Disc,MD0,Oiter);
        end
    end
end

if writedata
    fclose(fid);
end        

end

%{
oa = importdata('OA from web/MA.36.3.7.6.3.finney.txt');
oa = oa(:,5:10);
[N,n] = size(oa);
q = zeros(n,1);
for i = 1:n
    q(i) = length(unique(oa(:,i)));
end
[aveCD,A,BB,VEC] = aveCD_LevelPerm(oa,q);

InIter = 200;
OutIter = 200*InIter;

[MD0, MDvec, D0] = TA_uniform_FF_design_MD_abandon(oa,q,OutIter,InIter,1e-2,1e-5,10,0);
D0 = sortrows(D0);
fprintf('%.8f\n',MDvec);
%}