function [D0,E0,id0] = DWPsub_from_OA(D,q0,DiscMeasure,iswrite,filename)
%20150621 by Xiaopang
%���� Discrepancy induced pattern �����С��ԭ��, �� D ��ѡȡ�����D0,
%������Ƶĸ���ˮƽ��Ϊq0(i).
%������Ҫ��D�и���ˮƽ����������D0��������һ��
%INPUT:
%   D: һ���ǶԳ� OA
%   q0: n-by-1 ��������¼OA���е�ˮƽ��
%   DiscMeasure: 'CD','WD','MD'
%   iswrite: �Ƿ񽫹��̽��д���ļ��� 0��default����д��1��д
%OUTPUT:
%   D0: �õ������
%   E0: �õ�����Ƶ�pattern (E_1(D0),...,E_n(D0))
%   id0: D(:,id0) = D0;

epsilon = 1e-10;
if nargin < 4 || isempty(iswrite)
    iswrite = 0;
end
if nargin < 3
    DiscMeasure = 'CD';
end
if ~strcmp(DiscMeasure,'CD') && ~strcmp(DiscMeasure,'WD') && ~strcmp(DiscMeasure,'MD')
    error('DiscMeasure must be CD,WD or MD!\n');
end
n = size(D,2);
q = zeros(n,1);
for i = 1:n
    q(i) = length(unique(D(:,i)));
end
if ~issorted(q)
    [q,id] = sort(q);
    D = D(:,id);
end
uq = unique(q); %���г��ֵĲ�ͬˮƽ��
m = histc(q,[uq-0.5;inf]); 
m(end) = []; %��ˮƽ���ӵĸ�����m(i) Ϊˮƽ��Ϊ uq(i) �����Ӹ���
if ~issorted(q0)
    q0 = sort(q0);
end
uq0 = unique(q0); %�������ˮƽ�����
m0 = histc(q0,[uq0-0.5;inf]); 
m0(end) = []; %�������ˮƽ��Ϊ uq(i) �����Ӹ���Ϊ m(i)
if length(m0)~=length(m) || any(uq0~=uq) || any(m0>m)
    error('ÿ��ˮƽ�������Ӷ���Ҫ���֣�\n');
end

if iswrite
    %outfile = fopen('MA LHD 20151102 MD oa36 comp order/result.txt','w');
    if nargin < 5
        filename = 'result.txt';
    end
    outfile = fopen(filename,'w');
end

dd = zeros(1,m0(1));
for i = 2:length(m0)
    dd = cat(2,dd,(dd(end)+m(i-1))*ones(1,m0(i))); 
end
id = 1:m0(1);
for i = 2:length(m)
    id = cat(2,id,(1:m0(i)));
end

E0 = inf(1,length(q0));
while id(end)~=-1
    if ~iswrite
        [E,y] = Weighted_wordtype_pattern(D(:,id+dd),q0,DiscMeasure);
    else
        [E,y,A,BB] = Weighted_wordtype_pattern(D(:,id+dd),q0,DiscMeasure);
    end
    if iswrite
        fprintf(outfile,'%.8e ',y);
        fprintf(outfile,'%d ',id+dd); fprintf(outfile,'%.16e ',E); 
        fprintf(outfile,'%.16e ',A); fprintf(outfile,'%.16e ',BB); 
        fprintf(outfile,'\n');
        %{
        fprintf('%.4e ',y); fprintf('%d ',id+dd); fprintf('\n'); 
        fprintf('%.4e ',E); fprintf('\n');
        fprintf('%.4e ',A); fprintf('\n');
        fprintf('\n');
        %}
    end
    vv = find(abs(E-E0)>epsilon,1);
    if ~isempty(vv) && E(vv)<E0(vv)
        E0 = E;
        id0 = id;
    end
    id = choosenextgroup(m,m0,id); 
end
id0 = id0+dd;
D0 = D(:,id0);

if iswrite
    %{
    fprintf(outfile,'\n The best subarray is:');
    fprintf(outfile,'%d ',id0); fprintf(outfile,'\n');
    fprintf(outfile,'%.8f ',E0); fprintf(outfile,'\n');
    fprintf('\n The best subarray is:');
    fprintf('%d ',id0); fprintf('\n');
    fprintf('%.8f ',E0); fprintf('\n');
    %}
    fclose(outfile);
end

end









%{
oa = importdata('OA from web/MA.16.2.6.4.3.txt');
q0 = [2;2;4;4];
[D0,E0,id0] = DWPsub_from_OA(oa,q0,'WD',0);

oa = importdata('OA from web/MA.36.3.12.2.11.txt');
[N,ncols] = size(oa);
q = zeros(ncols,1);
for i = 1:ncols
    q(i) = length(unique(oa(:,i)));
end
[q,id] = sort(q);
oa = oa(:,id);
q0 = [2;2;3];
[D0,E0,id0] = DWPsub_from_OA(oa,q0,'WD',0);
%}

