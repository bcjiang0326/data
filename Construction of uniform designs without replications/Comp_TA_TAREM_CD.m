function Comp_TA_TAREM_CD(InIter,lambda)
%{
data = [ 3,3,18; 3,4,42; 3,5,48; 4,2,12; 4,3, 44; 
    3,5,102; 3,5,201; 3,6,201;3,6,502;5,4,100; 5,4, 400];
%}
%data = [3,5,102; 3,5,201; 3,6,102; 3,6,201;3,6,300; 3,6,402; 3, 6,501;5,4,100; 5,4,200; 5,4, 300; 5,4, 400; 5,4,500];
%
nzone = (102:102:2142)';
ln = length(nzone);
data = [ones(ln,1)*3,ones(ln,1)*7,nzone];
%}
%{
nzone = [36:33:714]';
ln = length(nzone);
data = [ones(ln,1)*3,ones(ln,1)*6,nzone];
%}

fid = fopen('TAREM_TA_CD.txt','a');
fprintf(fid,'InIter= %d, lambda = %d\n',InIter,lambda);
for i = 1:size(data,1)
    cout_d(data(i,:));
    levels = data(i,1); s = data(i,2); n = data(i,3);
    
    q = levels*ones(s,1);
    N = levels^s;
    
    % 公用的参数
    writedata = 0; Reps = 30;
    
    %modified TA
    T0 = 0.01;
    T1 = 1e-6;
    %lambda = 100;
    %InIter = 100;
    OutIter = InIter*lambda;
    t1 = cputime;
    if n > 0.5*N && (mod(levels,2)~=0 || s==2 || levels==2)
        [Disc1, Disc1vec,ID1] = TA_REM_CD_f(N-n,q,OutIter,InIter,T0,T1,Reps,writedata);
    else
        [Disc1, Disc1vec,ID1] = TA_REM_CD_f(n,q,OutIter,InIter,T0,T1,Reps,writedata);
    end
    t1 = cputime-t1;
    u1 = [Disc1, mean(Disc1vec), std(Disc1vec), t1];
    if n > 0.5*N && (mod(levels,2)~=0 || s==2 || levels==2)
        [~, ~,~,compDisc1vec] = TA_REM_CD_f(N-n,q,OutIter,InIter,T0,T1,Reps,writedata);
        u1 = [min(compDisc1vec), mean(compDisc1vec), std(compDisc1vec), t1];
    end
    cout_8f(u1);
    
    %TA
    t2 = cputime;
    [Disc2, Disc2vec,D2] = TA_CD(n,q,OutIter,InIter,T0,T1,Reps,writedata);
    t2 = cputime-t2;
    u2 = [Disc2, mean(Disc2vec), std(Disc2vec), t2];
     [~, ~, ~, K] = TA_CD(n,q,OutIter,InIter,T0,T1,Reps,writedata);
    cout_8f(u2);
    
    u0 = [i,n,s,levels,u1(2:end),u2(2:end),K];
     if n > 0.5*N && (mod(levels,2)~=0 || s==2 || levels==2)
         fprintf('%.6f(%.2e) & *%.1f & %.6f(%.2e) & %.1f & %.2f \n',u0(5:end));
         fprintf(fid,'%d & %d & %d  %d & %.6f(%.2e) & $\\ast$ %.1f & %.6f(%.2e) & %.1f & %.2f\\\\ \n',u0);
     else
         fprintf('%.6f(%.2e) & %.1f & %.6f(%.2e) & %.1f & %.2f \n',u0(5:end));
         fprintf(fid,'%d & %d & %d  %d & %.6f(%.2e) & %.1f & %.6f(%.2e) & %.1f & %.2f\\\\ \n',u0);
     end
end

fclose(fid);
end