%clear;clc;


Disc_fun = @CD2_value;

oa = importdata('OA from web/MA.36.3.12.2.11.txt');
q0 = [2;2;2;3;3;3];
[D0,A0,id0] = MA_from_asymOA(oa,q0,0);
[D1,E1,id1] = DWPsub_from_OA(oa,q0,'CD',0);
E0 = Disc_wordlength_pattern(D0,[],'CD');
A1 = GWP_asym(D1);

[N,n] = size(D0);
q = N*ones(n,1);
K = 1e4; 
vec0 = zeros(K,1); 
vec1 = zeros(K,1); 
for i = 1:K
    L1 = rand_oa_lhd(D0,q0,1);
    L2 = rand_oa_lhd(D1,q0,1);
    vec0(i) = Disc_fun(L1,q);
    vec1(i) = Disc_fun(L2,q);
end


ymin=min([vec0;vec1]);
ymax=max([vec0;vec1]);
x=linspace(ymin,ymax,100); %将最大最小区间分成20个等分点(19等分),然后分别计算各个区间的个数 
yy0=hist(vec0,x)/K;  %计算各个区间的个数
yy1=hist(vec1,x)/K;  %计算各个区间的个数



figure;
%plot(x,yy0,'b-',x,yy1,'b--'); hold on;
bar(x,[yy0;yy1]','grouped'); 
legend('D_0-based','D_1-based');
hold off;




