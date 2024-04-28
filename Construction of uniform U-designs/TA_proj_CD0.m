function [ptn0, ptn_vec, D0] = TA_proj_CD0(dim, N,q,OutIter,InIter,T0,T1,Reps,writedata)
% 20131001 ����汾��dim ֻ��Ϊ����
% ���� TA �㷨����ض�ά�ȣ�dim����ͶӰ������
% ��������ƣ��������� columnwise-pairwise ����
% ������Hickernell and Liu (2002) ��� projection discrepancy pattern,
%       ����������ر�ά��ͶӰ���������ŵ���ơ������ܲ��� TA �㷨��
%       .
% INPUT:
%       dim��scalar ͶӰά��
%       N: run-size
%       q: n-by-1 vector, ��Ƹ����ӵ�ˮƽ��
%       OutIter: �ⲿ��������
%       InIter: �ڲ���������
%       T0: ��ʼ��ֵ��λ��(0,1)֮�䣬����1e-2
%       T1: ������ֵ��λ��(0,1)֮�䣬����1e-6
%       Reps: �ظ����� (default,1)
%       writedata: �Ƿ񽫹�������д���ļ��� 1,д��; 0(default) ��д
% OUTPUT:
%       ptn0: �õ����ž���� Squared Discrepancy
%       ptn_vec: Reps-by-1 vector, ��¼ÿ���ظ������Ž�
%       D0: �õ������ž���

if min(size(q)) ~= 1 || N ~= floor(N) || any( N./q ~= floor(N./q) )
    error(' Input variable q must be a s-by-1 vector!\n');
end
if nargin < 8
    Reps = 1;
    writedata = 0;
elseif nargin < 9
    writedata = 0;
end
if Reps > 1
    writedata = 0;
end
n = length(q);

% �����ļ����
if writedata
    fid = fopen('T_Disc_Disc0.txt','w');
end

% ���ù������
epsilon = 1e-12;

% ���ȼ�����ֵ�仯�ı���,��ÿ InIter �ε�����
% ����һ����ֵ T = T*ratio     
ratio = (T1/T0)^(1/(ceil(OutIter/InIter)-1));


ptn0 = inf; 
ptn_vec = zeros(Reps,1);

for rep = 1:Reps
    % ��������һ��������ظ��ľ���
    D = rand_U_Type_design(q,N);
    DD = (D+0.5)./(2*ones(N,1)*q')-0.25;
    aD = abs(DD);
    
    % ��ʼ�� PDisc
    ptn = CD2_pattern(D,q,dim);

    % �жϵ����ظ��ĳ�ʼ���Ƿ����
    ptn_vec(rep) = ptn;
    if ptn < ptn0-epsilon
        ptn0 = ptn;
        if nargout > 2
            D0=D;
        end
    end
    
    %�����㷨
    T = T0;
    for Oiter = 1:OutIter
        if mod(Oiter,InIter)==0
            T = T*ratio;
        end
        
        %����Ҫ�������м���������Ԫ��    
        k0 = randi( n ); % ���ѡ�� k0 ��
        i0 = randi( N );
        j0 = randi( N );
        while j0==i0 || D(i0,k0) == D(j0,k0)    
            j0=randi( N );
        end
        
        %i0 = 2; j0 = 4; k0 = 1;
        
   
        %����delta(kk) = delta1 + delta21 + delta22
        delta1 = 0; delta21 = 0; delta22 = zeros(1,N);
        id  = 1:dim-1;
        while id(1)~=-1
            id0 = id;
            id0(id0>=k0) = id0(id0>=k0)+1;
            delta1 = delta1+prod(aD(i0,id0)-2*DD(i0,id0).^2)-prod(aD(j0,id0)-2*DD(j0,id0).^2);
            delta21 = delta21+prod(aD(i0,id0))-prod(aD(j0,id0));
            for t = 1:N
                if t~=i0 && t~=j0
                    delta22(t) = delta22(t)+prod( aD(i0,id0)+aD(t,id0)-abs(DD(i0,id0)-DD(t,id0)) ) ...
                        -prod( aD(j0,id0)+aD(t,id0)-abs(DD(j0,id0)-DD(t,id0)) );
                end
            end
            id = nchoosek_next(id,n-1,dim-1);
        end
        delta1 = delta1*(-2/N)*(aD(j0,k0)-aD(i0,k0))*(1-2*(aD(j0,k0)+aD(i0,k0)));
        delta21 = delta21*2^dim/N^2*(aD(j0,k0)-aD(i0,k0));
        delta22 = delta22*(aD(j0,k0)-aD(i0,k0)-abs(DD(j0,k0)-DD(:,k0))+abs(DD(i0,k0)-DD(:,k0)))*2/N^2;
        delta = delta1+delta21+delta22;
        
        %�ж��Ƿ����
        %if delta <  PDisc*T-epsilon
        if  (rand()<0.001 && T < 1e-4) ||  delta<  ptn*T-epsilon
            % ���� D,DD,PDisc
            temp=D(i0,k0);D(i0,k0)=D(j0,k0);D(j0,k0)=temp;
            temp=DD(i0,k0);DD(i0,k0)=DD(j0,k0);DD(j0,k0)=temp;
            temp=aD(i0,k0); aD(i0,k0)=aD(j0,k0); aD(j0,k0)=temp;
            ptn = ptn + delta;
            % �� delta<0 ˵����ʼ�½�.
            if delta < -epsilon
                % �������ȱ��θ��Ž⣬���¼֮
                if ptn < ptn_vec(rep)-epsilon
                    ptn_vec(rep) = ptn;
                    % ������ȫ�ָ��Ž⣬���¼֮
                    if ptn_vec(rep) < ptn0
                        ptn0 = ptn_vec(rep);
                        if nargout > 2
                            D0 = D;
                        end
                    end
                end
            end
        end
        if writedata
            fprintf(fid,'%.8f %.8f %.8f %.8f %.8f\n',delta/ptn,T,ptn,ptn0,Oiter);
        end
    end
end

if writedata
    fclose(fid);
end
end

%{
% ���Ժ���
D=[2 4 2 4; 4 3 2 1; 3 4 4 2; 2 1 3 1; 4 2 3 4; 3 1 1 3; 1 2 1 2; 1 3 4 3];
D = sortrows(D-1); % Uniform (8,4^4)-design

s = max(max(D))+1;
[N,n] = size(D);
q = s*ones(n,1);


InIter = 200;
OutIter = 200*InIter;
Reps = 1;
dims = [2,3];
[ptn0, ptn_vec, D0] = TA_proj_CD(dims,N,q,OutIter,InIter,1e-2,1e-6,Reps,0);
fprintf('%.4e ',ptn0);fprintf('\n');
%}