function Comp_SAIPM_SAREM(isWD)
%20130124
%INPUT:
%   isWD: 0,CD; 1,WD.
%{
data = [3 5 48; 3 5 102; 3 5 201; 3 7 48; 3 7 201; 3 7 1002;
    4 4 36; 4 4 100; 4 4 200; 4 5 36; 4 5 200; 4 5 500;
    5 4 100; 5 4 200; 5 4 500; 5 5 200; 5 5 500; 5 5 1000; 
    6 4 102; 6 4 204; 6 4 504; 3 8 201; 3 8 501; 3 8 1002;
    4 6 200; 4 6 500; 4 6 1000; 3 10 201; 3 10 501; 3 10 1002];
%}
data = [5 4 50];

Reps = 30;

if isWD
    fname = strcat('SAIPM_SAREM_WD0226','.txt');
else
    fname = strcat('SAIPM_SAREM_CD0226','.txt');
end
fid = fopen(fname,'w');
for cs = 1:size(data,1)
    cout_d(data(cs,:));
    levels = data(cs,1); s = data(cs,2); n = data(cs,3);
    
    q = levels*ones(s,1);
    N = levels^s;

if cs < 27
    if isWD
        T0 = 1/levels; % WD 情形
        beta_in = 0.99;
    else
        T0 = 1/(0.5*levels); %CD 情形
        beta_in = 0.99;
    end
    
    InIter = 10;
    if N < 500
        OutIter = 10;
        beta_out = 0.9;
    else
        OutIter = 2;
        beta_out = 0.8;
    end
    t1 = cputime;
    [Disc1,Disc_vec1] = SA_IPM(n,q,OutIter,InIter,T0,beta_in,beta_out,isWD,Reps);
    t1 = cputime-t1;
    cout_8f([mean(Disc_vec1),std(Disc_vec1),t1]);
end
    
    if isWD
        T0 = 1e-4;
        T1 = 1e-7;
    else
        T0 = 1e-4;
        T1 = 1e-7;
    end
    beta_in = 0.8;
    InIter = 10;
    OutIter = 500;
    
    if isWD
        %{
        if n > 0.5*N %&& (mod(levels,2)==1 || s == 2)
            useComp = 1;
            t2 = cputime;
            [Disc2,Disc_vec2] = SA_REM_WD(N-n,q,OutIter,InIter,T0,T1,beta_in,Reps);
            Disc2 = comp_WD_value(Disc2,q,N-n);
            t2 = cputime-t2;
        else
        %}
            %useComp = 0;
            t2 = cputime;
            [Disc2,Disc_vec2] = SA_REM_WD(n,q,OutIter,InIter,T0,T1,beta_in,Reps);
            t2 = cputime-t2;
        %end
    else
        %{
        if n > 0.5*N && (mod(levels,2)==1 || s == 2)
            useComp = 1;
            t2 = cputime;
            [Disc2,Disc_vec2] = SA_REM_CD(N-n,q,OutIter,InIter,T0,T1,beta_in,Reps);
            Disc2 = comp_WD_value(Disc2,q,N-n);
            t2 = cputime-t2;
        else
        %}
            %useComp = 0;
            t2 = cputime;
            [Disc2,Disc_vec2] = SA_REM_CD(n,q,OutIter,InIter,T0,T1,beta_in,Reps);
            t2 = cputime-t2;
        %end
    end
    
    
    cout_8f([mean(Disc_vec2),std(Disc_vec2),t2]);
    if cs < 28
        u = [cs,n,s,levels,N,mean(Disc_vec2),std(Disc_vec2),t2,mean(Disc_vec1),std(Disc_vec1),t1];
        fprintf(fid,'%d & %d & %d %d & %d & %.6f & %.2e & %.2f & %.6f & %.2e & %.2f \\\\ \n',u);
    else
        u = [cs,n,s,levels,N,mean(Disc_vec2),std(Disc_vec2),t2];
        fprintf(fid,'%d & %d & %d %d & %d & %.6f & %.2e & %.2f &     &     &    \\\\ \n',u);
    end
end

fclose(fid);
end