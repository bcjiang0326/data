%对于 MA.72.4.1.6.5.3.8.2.7 形为 MA.72.4.1.6.n1.3.n2. 的 MAD, GMA 子表不同的情况，
%本程序搜索这些子表在水平置换意义下的左右均匀设计，水平置换过程采取随机算法。
oa = importdata('../Web OA/MA.72.4.1.6.5.3.8.2.7.txt');
n1n2=[1,2;2,2;2,3;2,4;2,6;3,2;3,3;3,4;3,5;3,6;3,7;3,8;4,1;4,3;4,4;
    4,5;4,6;5,3;5,4;5,5;5,6;5,7];

flag = 'WD';
Reps = 10;
InIter = 100;
AllIter = 100*InIter;
for i = 1:size(n1n2,1)
    %n1 = n1n2(i,1); n2 = n1n2(i,2);
    n1 = 4; n2 = 1;
    q0 = [4;6*ones(n1,1);3*ones(n2,1)];
    T1 = 0; T2 = 0;
    
    %阶段一: 构造 MAD 和 GMA 设计
    t = cputime;
    [~,id1,aveDisc1] = MADsub_from_OA(oa,q0,flag);
    T1 = T1+cputime-t;
    t = cputime;
    [~,id2,A2] = GMAsub_from_OA(oa,q0);
    T2 = T2+cputime-t;
    A1 = GWP_asym(oa(:,id1),q0);
    aveDisc2 = aveDisc_LevelPerm(oa(:,id2),q0,flag);
    if A1(3)<A2(3)
        fprintf('Warning!\n');
    end
    
    %阶段二：实施水平置换
    t = cputime;
    [WD1, WD1vec, D1] = TA_uniform_FF_design_WD(oa(:,id1),AllIter,InIter,1e-3,1e-6,Reps);
    T1 = T1+cputime-t;
    t = cputime;
    [WD2, WD2vec, D2] = TA_uniform_FF_design_WD(oa(:,id2),AllIter,InIter,1e-3,1e-6,Reps);
    T2 = T2+cputime-t;
    
    
    y = [n1,n2,WD1,std(WD1vec),A1(3),T1,WD2,std(WD2vec),A2(3),T2];
    fprintf('%d & %d & %.6f & %.1e & %.1f & %.2f && %.6f & %.1e & %.1f & %.2f\n',y);
end