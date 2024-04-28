function Comp_SAIPM_TAREM( )

%data = [3 10 201; 3 10 501; 3 10 1002];
%
nzone = [102:102:2142]';
ln = length(nzone);
data = [ones(ln,1)*3,ones(ln,1)*7,nzone];
%}
%{
nzone = [102;2040;2142];
ln = length(nzone);
data = [ones(ln,1)*3,ones(ln,1)*7,nzone];
%}

%data = [3,5,201; 3,6,201; 3,7,201; 3,8,201];
%fid = fopen('TAREM_SAIPM_WD.txt','w');
for i = 1:size(data,1)
    cout_d(data(i,:));
    levels = data(i,1); s = data(i,2); n = data(i,3);
    
    q = levels*ones(s,1);
    N = levels^s;
    
    % 公用的参数
    writedata = 0; Reps = 30;
    
    %{
    %modified TA
    T0 = 0.01;
    T1 = 1e-6;
    %{
    if n > 0.5*N
        InIter = min(floor((N-n)/2),100);
    else
         InIter = min(floor(n/2),100);
    end
    %}
    InIter = 100;
    lambda = 100;
    OutIter = InIter*lambda;
    t1 = cputime;
    if n > 0.5*N
        [Disc1, Disc1vec,ID1] = TA_REM_WD_f(N-n,q,OutIter,InIter,T0,T1,Reps,writedata);
    else
        [Disc1, Disc1vec,ID1] = TA_REM_WD_f(n,q,OutIter,InIter,T0,T1,Reps,writedata);
    end
    t1 = cputime-t1;
    u1 = [Disc1, mean(Disc1vec), std(Disc1vec), t1];
    if n > 0.5*N
        u1(1) = WD2_From_comp(q,n,u1(1));
        u1(2) = WD2_From_comp(q,n,u1(2));
        u1(3) = (N-n)^2/n^2*u1(3);
    end
    cout_8f(u1);
    %}
    
    %{
    %TA
    InIter = min(floor(n/2),100);
    OutIter = InIter*lambda;
    t2 = cputime;
    [Disc2, Disc2vec,D2] = TA_WD(n,q,OutIter,InIter,T0,T1,Reps,writedata);
    t2 = cputime-t2;
    u2 = [Disc2, mean(Disc2vec), std(Disc2vec), t2];
    cout_8f(u2);
    %}
    
    %
    % SA-IPM
    Tini = 1/levels;
    Tf = 0.99;
    InIter = 10;
    if N < 500
        OutIter = 10;
        Tr = 0.9;
    else
        OutIter = 2;
        Tr = 0.8;
    end
    t3 = cputime;
    [Disc3, Disc3vec, ID3] = SA_IPM(n,q,OutIter,InIter,Tini,Tf,Tr,1,Reps,writedata);
    t3 = cputime-t3;
    u3 = [Disc3, mean(Disc3vec), std(Disc3vec), t3];
    cout_8f(u3);
    %}
    
    %{
    u = [i,n,s,levels,u1(2:end),u3(2:end)];
    if n<0.5*N
        fprintf(fid,'%d & %d & %d  %d & %.6f(%.2e) & %.1f~ & %.6f(%.2e) & %.1f \\\\ \n',u);
    else
        fprintf(fid,'%d & %d & %d  %d & %.6f(%.2e) & %.1f$^{\\ast}$ & %.6f(%.2e) & %.1f \\\\ \n',u);
    end
    %}
    %{
        u = [i,n,s,levels,u1(2:end),u2(2:end),u3(2:end)];
        fprintf('%.6f(%.2e) & %.1f & %.6f(%.2e) & %.1f & %.6f(%.2e) & %.1f\\\\ \n',u(5:end));
        fprintf(fid,'%d & %d & %d  %d & %.6f(%.2e) & %.1f & %.6f(%.2e) & %.1f & %.6f(%.2e) & %.1f\\\\ \n',u);
     %}
end
% fclose(fid);
end