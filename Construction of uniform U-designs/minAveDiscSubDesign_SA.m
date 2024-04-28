function [aveDisc0,cols0,WB0,A0,D0] = minAveDiscSubDesign_SA(D,q0,DiscMeasure,OutIter,InIter,T0,T1,writedata)
% 20160518 Xiaopang
% ���� TA �㷨��ȡ D �и���ˮƽΪ q0 ������ƣ�ʹ�����ɵ�����ƾ��н�С��ƽ��������
%INPUT:
%   D: һ���ǶԳ� OA
%   q: n-by-1 ��������¼OA���е�ˮƽ��
%   DiscMeasure: 'CD','WD','MD'
%   OutIter: �ⲿ��������
%   InIter: �ڲ���������
%   T0: ��ʼ��ֵ��λ��(0,1)֮�䣬����1e-2
%   T1: ������ֵ��λ��(0,1)֮�䣬����1e-5
%   writedata: �Ƿ񽫹�������д���ļ� -- 1,д��; 0(default)
%OUTPUT:
%   aveDisc: ƽ��������
%   cols: ����Ƶ��б��
%   WB: discrepancy induced weighted wordtype pattern
%   A: generalized wordlength pattern
%   D0: �����
n = size(D,2);
if nargin < 3
    DiscMeasure = 'CD';
end
if nargin < 4
    OutIter = 10^3;
end
if nargin < 5
    InIter = 100;
end
if nargin < 6
    T0 = 1e-2;
end
if nargin < 7
    T1 = 1e-6;
end
if nargin < 8
    writedata = 0;
end

% ���ȼ�����ֵ�仯�ı���,��ÿ InIter �ε�����
% ����һ����ֵ T = T*ratio
if T0 ~= 0
    ratio = (T1/T0)^(1/(ceil(OutIter/InIter)-1));
else
    ratio = 1;
end
epsilon = 1e-11;

q = zeros(n,1);
for i = 1:n
    q(i) = length(unique(D(:,i)));
end
if size(q,2)~=1 || size(q,1)~=n
    error('q must be an n-by-1 vector!\n');
end
if ~issorted(q)
    [q,id] = sort(q);
    D = D(:,id);
end
if ~strcmp(DiscMeasure,'CD') && ~strcmp(DiscMeasure,'WD') && ~strcmp(DiscMeasure,'MD')
    error('DiscMeasure must be CD,WD or MD!\n');
end
uq = unique(q); %���г��ֵĲ�ͬˮƽ��
m = histc(q,[uq-0.5;inf]);
m(end) = []; %��ˮƽ���ӵĸ�����m(i) Ϊˮƽ��Ϊ uq(i) �����Ӹ���

q0 = sort(q0);
uq0 = unique(q0); %������г��ֵĲ�ͬˮƽ��
m0 = histc(q0,[uq0-0.5;inf]);
m0(end) = [];

% �������
if writedata
    fid = fopen('aveDisc.txt','w');
end


if ~all(ismember(uq0,uq))
    error('The elements of q0 must be the number of certain factor in D !\n');
end

for i = 1:length(uq)
    if ~ismember(uq(i),uq0)
        m0 = [m0(1:i-1),0,m0(i:end)];
    end
end

for i = 1:length(uq)
    if  m0(i) > m(i)
        error('Wrong q0!\n');
    end
end

%�������һ����ʼ�������
cum = cat(1,0,cumsum(m));
cols = zeros(1,0);
for i = 1:length(m0)
    if m0(i)==0
        continue;
    end
    subcols = randperm(m(i),m0(i));
    cols = sort(cat(2,cols,subcols+cum(i)));
end
%cols = [1,2,4,17,18,19];
[WB,aveDisc,A] = Weighted_wordtype_pattern(D(:,cols),q0,DiscMeasure);
if length(cols) == n
    return;
end

aveDisc0 = aveDisc;
WB0 = WB;
A0 = A;
cols0 = cols;

%�� groups_in �� groups_out ��¼ÿ��ˮƽ����������е��кͲ���������е��С�
cum0 = cat(1,0,cumsum(m0));
groups_in = cell(1,length(uq));
groups_out = cell(1,length(uq));
for k = 1:length(uq)
    bg = cum(k)+1;
    if k == length(uq)
        ed = n;
    else
        ed = cum(k+1);
    end
    bg0 = cum0(k)+1;
    if k == length(uq)
        ed0 = length(q0);
    else
        ed0 = cum0(k+1);
    end
    groups_in{k} = cols(bg0:ed0);
    groups_out{k} = setdiff(bg:ed,groups_in{k});
end

Tini = 1/q(end);
beta_in = 0.99;
beta_out = 0.9;
for Oiter = 1:OutIter
    T=Tini;
    ct=0;
    while ct < InIter
        ct = ct + 1;
        
        %�������һ���µ������
        k = randi(length(uq)); %����ѡ��һ����
        while m0(k) == 0 || m0(k) == m(k)
            k = randi(length(uq));
        end
        i = randi(m0(k)); %groups_in(k) �ĵ�i��Ԫ�س���
        j = randi(m(k)-m0(k));%groups_out(k) �ĵ�j��Ԫ������
        cols_ = cols;
        cols_(cum0(k)+i) = groups_out{k}(j);
        %��������Ƶ�ƽ��������
        [WB_,aveDisc_,A_] = Weighted_wordtype_pattern(D(:,cols_),q0,DiscMeasure);
        
        delta = aveDisc_-aveDisc;
        %�ж��Ƿ����
        if delta < -epsilon || rand() < exp(-delta/aveDisc/T)
            WB = WB_;
            aveDisc = aveDisc_;
            A = A_;
            cols = cols_;
            kk = groups_in{k}(i);
            groups_in{k}(i) = groups_out{k}(j);
            groups_out{k}(j) = kk;
            if aveDisc < aveDisc0-epsilon
                aveDisc0 = aveDisc;
                WB0 = WB;
                A0 = A;
                cols0 = cols;
            end
        end
        if writedata
            fprintf(fid,'%.8f\n',aveDisc);
            fprintf(fid,'%.3e ',WB);
            fprintf(fid,'\n');
            fprintf(fid,'%.3e ',A);
            fprintf(fid,'\n');
        end
        if ct == 1
            fprintf('%d, %.8f %.2f\n',Oiter,aveDisc0,min(1,exp(-delta/aveDisc/T)))
        end
        T = T*beta_in;
    end
    Tini = Tini*beta_out;
end
if writedata
    fclose(fid);
end
if nargout > 4
    D0 = D(:,cols);
end
end
