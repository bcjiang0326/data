%本程序计算阶段1选取MAD或GMA子表所用时间
oa = importdata('../Web OA/MA.36.2.11.3.12.txt');
n1n2_ = [3,3; 5,2; 5,3; 5,5; 4,6; 9,3; 3,8; 3,9; 3,10; 5,9];
n1n2 = zeros(size(n1n2_));
n1n2(:,1) = n1n2_(:,2); n1n2(:,2) = n1n2_(:,1);
allDisc = {'WD';'CD';'MD'};
fid = fopen('Stage1Time_MA.36.3.12.2.11.txt','w');
for k = 1:size(n1n2,1)
    n1 = n1n2(k,1); n2 = n1n2(k,2);
    q0 = [2*ones(n1,1);3*ones(n2,1)];
    
    T1 = zeros(1,3);
    for j = 1:3
        flag = allDisc{j};
        T1(j) = cputime;
        [D0,id0,aveDisc0] = MADsub_from_OA(oa,q0,flag);
        T1(j) = cputime-T1(j);
        fprintf('(%d,%d):',n1,n2);
        fprintf('%d ',id0);
        fprintf('& ');
    end
    
    T2 = cputime;
    [D1,id1,A1] = GMAsub_from_OA(oa,q0);
    T2 = cputime-T2;
    fprintf(fid,'Stage2Time(%d,%d): %.2f  %.2f  %.2f  %.2f\n',n1,n2,T1,T2);    
    fprintf('%d ',id1);
    fprintf('\n');
end
fclose(fid);