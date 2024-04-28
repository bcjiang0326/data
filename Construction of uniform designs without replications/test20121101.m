
%INPUT:
%   isWD: 0,CD; 1,WD.
data = [3 10 201; 3 10 501; 3 10 1002];
%}
%data = [3 3 15];

isWD = 1;
Reps = 30;

if isWD
    fname = strcat('SAIPM_SAREM_WD0222','.txt');
else
    fname = strcat('SAIPM_SAREM_CD0222','.txt');
end
fid = fopen(fname,'a');
for cs = 1:size(data,1)
    cout_d(data(cs,:));
    levels = data(cs,1); s = data(cs,2); n = data(cs,3);
    
    q = levels*ones(s,1);
    N = levels^s;
    
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
        t2 = cputime;
        [Disc2,Disc_vec2] = SA_REM_WD(n,q,OutIter,InIter,T0,T1,beta_in,Reps);
        t2 = cputime-t2;
    else
        t2 = cputime;
        [Disc2,Disc_vec2] = SA_REM_CD(n,q,OutIter,InIter,T0,T1,beta_in,Reps);
        t2 = cputime-t2;
    end
    
    
    cout_8f([Disc2,mean(Disc_vec2),t2]);
    
    u = [cs,n,s,levels,N,mean(Disc_vec2),std(Disc_vec2),t2];
    fprintf(fid,'%d & %d & %d %d & %d & %.6f & %.2e & %.2f &     &     &    \\\\ \n',u);
end

fclose(fid);
