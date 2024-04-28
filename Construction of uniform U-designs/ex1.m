%{
D = importdata('OA from web/MA.36.3.12.2.11.txt');
[N,n] = size(D);
q = zeros(n,1);
for i = 1:n
    q(i) = length(unique(D(:,i)));
end
[q,id] = sort(q);
D = D(:,id);
q0 = [2;2;2;3;3;3];

results = importdata('result.txt');
Y = results(:,1);
cols = results(:,2:7);
Ep = results(:,10:13);
clear('results');
A =  zeros(0,4);
BB = zeros(0,16);
for i = 1:size(cols,1)
    [gwp,BB_] = GWP_asym(D(:,cols(i,:)),q0);
    A = cat(1,A,gwp(3:6));
    BB = cat(1,BB,BB_');
end

m = [3;3];
vec = zeros(0,2);
for k = 1:size(BB,2)
    vec = cat(1,vec,Id2Design(m+1,k-1));
end
typeLen = sum(vec,2);
[typeLen,v] = sort(typeLen);
vec = vec(v,:);
BB = BB(:,v);
v = typeLen>2;
BB = BB(:,v);
vec = vec(v,:);
typeLen = typeLen(v);
%}
clear; clc;
load ex1.mat
k = 12; %保留到小数点后第k位
A = round(A*10^k)/10^k;
BB = round(BB*10^k)/10^k;
Ep = round(Ep*10^k)/10^k;
Y = round(Y*10^k)/10^k;

[BB,ib] = unique(BB,'rows'); % wordtype pattern (j >= 3 部分）
cols = cols(ib,:);
A = A(ib,:); % GWP
Ep = Ep(ib,:); % weighted wordtype pattern
Y = Y(ib,:); %averaged discrepancy


[Y,iy] = sort(Y);
cols = cols(iy,:);
BB = BB(iy,:);
A = A(iy,:);
Ep = Ep(iy,:);

[~,ia] = sortrows(A);
[~,iep] = sortrows(Ep);

%{
Nrep = 20;
for i = 1:Nrep
   x = [i,Y(i),iep(i),ia(i),cols(i,:)]; 
   %fprintf('%d & %.6f & %d & %d & %d %d %d %d %d %d \\\\\n',x);
   fprintf('%d & %.6f & %d & %d & %d %d %d %d %d %d &',x);
   fprintf('%.2f,',BB(i,:)); fprintf('\\\\\n');
end
%}

Disc_fun = @CD2_value;
k1 = 1; k2 = 600; 
D1 = D(:,cols(k1,:));
D2 = D(:,cols(k2,:));

[N,n] = size(D1);
q_lh = N*ones(n,1);
K = 1e5; 
vec1 = zeros(K,1); 
vec2 = zeros(K,1);
for i = 1:K
    L1 = rand_oa_lhd(D1,q0,1);
    vec1(i) = Disc_fun(L1,q_lh);
    L2 = rand_oa_lhd(D2,q0,1);
    vec2(i) = Disc_fun(L2,q_lh);
end

ymin=min([vec1;vec2]);
ymax=max([vec1;vec2]);
x=linspace(ymin,ymax,50); %将最大最小区间分成20个等分点(19等分),然后分别计算各个区间的个数 
yy1=hist(vec1,x);  %计算各个区间的个数
yy2=hist(vec2,x);  %计算各个区间的个数

figure;
plot(x,yy1,'k-o',x,yy2,'k-*'); hold on;
xlim([ymin,ymax]);
bar(x,[yy1;yy2]','grouped','k'); 
legend('D_1-based','D_2-based');
hold off;


