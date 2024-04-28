clear;clc;

N = 16; n = 5; s1 = 2; s2 = 4; s3 = 8;
Disc_fun = @CD2_value;
Fun = @CD_GWPoa_lhd;

oa16_8_2_3 = importdata('OA from web/oa.16.8.2.3.txt');
oa1 = oa16_8_2_3(:,1:n);
oa2 = importdata('OA from web/oa.16.5.4.2.a.txt');
oa3 = zeros(N,n);
for j = 1:n
    for lv = 0:s2-1
        v = find(oa2(:,j)==lv);
        oa3(v(1:2),j) = 2*lv;
        oa3(v(3:4),j) = 2*lv+1;
    end
end

[h,b,H,B] = Fun(N,n,s1);
A1 = GWP(oa1,s1);
Disc1 = h*[1,A1]'+b;
PD1 = H*[1,A1]'+B;

[h,b,H,B] = Fun(N,n,s2);
A2 = GWP(oa2,s2);
Disc2 = h*[1,A2]'+b;
PD2 = H*[1,A2]'+B;

[h,b,H,B] = Fun(N,n,s3);
A3 = GWP(oa3,s3);
Disc3 = h*[1,A3]'+b;
PD3 = H*[1,A3]'+B;

K = 1e5; 
q = N*ones(n,1);
vec1 = zeros(K,1); q1 = s1*ones(n,1);
vec2 = zeros(K,1); q2 = s2*ones(n,1);
vec3 = zeros(K,1); q3 = s3*ones(n,1);
for i = 1:K
    L1 = rand_oa_lhd(oa1,q1,1);
    L2 = rand_oa_lhd(oa2,q2,1);
    L3 = rand_oa_lhd(oa3,q3,1);
    vec1(i) = Disc_fun(L1,q);
    vec2(i) = Disc_fun(L2,q);
    vec3(i) = Disc_fun(L3,q);
end


ymin=min([vec1;vec2;vec3]);
ymax=max([vec1;vec2;vec3]);
x=linspace(ymin,ymax,100); %将最大最小区间分成20个等分点(19等分),然后分别计算各个区间的个数

figure;
subplot(1,3,1);  
yy1=hist(vec1,x)/K;  %计算各个区间的个数
bar(x,yy1) %画出概率密度分布图
ylim([0,0.4]);

subplot(1,3,2);
yy2=hist(vec2,x)/K;  %计算各个区间的个数
bar(x,yy2) %画出概率密度分布图
ylim([0,0.4]);

subplot(1,3,3);
yy3=hist(vec3,x)/K;  %计算各个区间的个数
bar(x,yy3) %画出概率密度分布图
ylim([0,0.4]);

figure;
bar(x,[yy1;yy2;yy3]','grouped');
legend('D_1-based','D_2-based','D_3-based');


figure;
v = x<0.048;
plot(x(v),yy1(v)*K,'b-',x(v),yy2(v)*K,'b--',x(v),yy3(v)*K,'b-*'); hold on;
legend('D_1-based','D_2-based','D_3-based');
bar(x(v),[yy1(v);yy2(v);yy3(v)]'*K,'grouped','b');
hold off;