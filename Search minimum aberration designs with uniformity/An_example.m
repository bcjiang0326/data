N = 32; n = 30; s = 4;
%N = 50; n = 48; s = 5;
Reps = 10000;
m = n/2; q = ones(m,1)*s;

fprintf('N = %d, n = %d, s = %d :',N,n,s);

file1 = strcat('MAD/WD/level',int2str(s),'/Permuted_MAD.',int2str(N),'.',int2str(s),'.',int2str(n),'.txt');
UD1 = importdata(file1);
mu1 = zeros(Reps,1);
for rep = 1:Reps
    %Step 1.先对D1的列进行随机置换
    UD1 = UD1(:,randperm(n));
    %Step 2.对每一列进行随机轮换
    v1 = ones(N,1)*randi(s,1,n);
    UD1 = mod(UD1+v1,s);
    %Step 3.生成设计 D1
    Da1 = (UD1(:,1:m)+1)/s;
    [~,Db1] = rand_oa_lhd(UD1(:,m+1:end),q,0);
    D1 = [Da1,Db1];
    %Step 4.抽样
    Y1 = Robot_arm_Sampling(D1);
    %Step 5.计算均值mu
    mu1(rep) = mean(Y1);
end
fprintf('%.6f  ',mean(mu1));

file2 = strcat('My UD/WD/level',int2str(s),'/UD.',int2str(N),'.',int2str(s),'.',int2str(n),'.txt');
UD2 = importdata(file2);
mu2 = zeros(Reps,1);
for rep = 1:Reps
    %Step 1.先对D1的列进行随机置换
    UD2 = UD2(:,randperm(n));
    %Step 2.对每一列进行随机轮换
    v2 = ones(N,1)*randi(s,1,n);
    UD2 = mod(UD2+v2,s);
    %Step 3.生成设计 D2
    Da2 = (UD2(:,1:m)+1)/s;
    [~,Db2] = rand_oa_lhd(UD2(:,m+1:end),q,0);
    D2 = [Da2,Db2];
    %Step 4.抽样
    Y2 = Robot_arm_Sampling(D2);
    %Step 5.计算均值mu2
    mu2(rep) = mean(Y2);
end
fprintf('%.6f\n',mean(mu2));