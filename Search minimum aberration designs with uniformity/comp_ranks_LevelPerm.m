%function [rate_a,Num] = comp_ranks_LevelPerm(n1,n2)
%20181128 By Bochuan Jiang
%�Ƚ��ӱ���GWP�µ��������ƽ���������µ�����
n1 = 5; n2 = 2;
D = importdata('OA from web/MA.36.3.7.6.3.finney.txt');
q = zeros(size(D,2),1);
for i = 1:size(D,2)
    q(i) = length(unique(D(:,i)));
end
if ~issorted(q)
    error('q is not sorted!\n');
end

uq = unique(q); %���г��ֵĲ�ͬˮƽ��
m = histc(q,[uq-0.5;inf]); 
m(end) = []; %��ˮƽ���ӵĸ�����m(i) Ϊˮƽ��Ϊ uq(i) �����Ӹ���
if length(uq) < 2
    error('Mat is not mixed!\n');
end
if n1>m(1) || n2 > m(2)
    error('Error!\n');
end
    
cols = enumerate_groups([m(1);m(2)],[n1;n2]);
cols(:,n1+1:end) = cols(:,n1+1:end)+m(1);



q0 = [3*ones(n1,1);6*ones(n2,1)];

N = size(cols,1); n = n1+n2;

fid = fopen('out.txt','w');
for i = 1:N
    [aveCD,A,BB] = aveCD_LevelPerm(D(:,cols(i,:)),q0);
    fprintf(fid,'%.12e ',A);
    fprintf(fid,'%.12e ',aveCD);
    fprintf(fid,'%.12e ',BB);
    fprintf(fid,'\n');
end
fclose(fid);

data = importdata('out.txt');
A = data(:,1:n);
aveCD = data(:,n+1);
BB = data(:,n+2:end);

%���þ���
k = 10; %������С������kλ
A = round(A*10^k)/10^k;
BB = round(BB*10^k)/10^k;
aveCD = round(aveCD*10^k)/10^k;

%�����ȡ����generalized wordtype pattern ��ͬ��subdesign
%�������Ӧ��minCD �Ƿ���ͬ
%��������ͬ��
%{
i = randi(length(aveDisc));
I = cols(i,:);
for j = 1:length(aveDisc)
    if j == i
        continue;
    end
    if norm(BB(j,:)-BB(i,:)) < 1e-8
        I = cat(1,I,cols(j,:));
    end
end
i1 = randi(size(I,1));
i2 = randi(size(I,1));
while i2==i1
    i2 = randi(size(I,1));
end
[CDmin1,CDvec1] = minCD_allPerms(D(:,cols(i1,:)));
[CDmin2,CDvec2] = minCD_allPerms(D(:,cols(i2,:)));
%}

%ȥ����������ͬ����
[BB,ib] = unique(BB,'rows');
cols = cols(ib,:);
A = A(ib,:);
aveCD = aveCD(ib,:);

%����ƽ������������
[aveCD,id] = sort(aveCD);
cols = cols(id,:);
A = A(id,:);
BB = BB(id,:);
[~,ia] = sortrows(A);
[~,ia] = sort(ia);

CDmin = zeros(size(cols,1),1);
for i = 1:size(cols,1)
    fprintf('i=%d\n',i);
    [CDmin(i),CDvec] = minCD_allPerms(D(:,cols(i,:)));
    if abs(aveCD(i)-mean(CDvec))>1e-8
        error('Wrong!\n');
    end
end
[~,iminCD] = sort(CDmin);
[~,iminCD] = sort(iminCD);

fid = fopen('rank.txt','w');
for i = 1:size(cols,1)
    fprintf(fid,'%d %d %d %d %d & ',cols(i,:));
    fprintf(fid,'%.8f & %.8f & ',CDmin(i),aveCD(i));
    fprintf(fid,'%d & %d %d \\\\ \n',i,ia(i),iminCD(i));
end
fclose(fid);

