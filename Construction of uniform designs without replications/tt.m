    levels = 27; s = 5; n = 27;
    q = levels*ones(s,1);
    
    % 公用的参数
    writedata = 0; Reps = 1;
    
    %modified TA
    T0 = 0.01;
    T1 = 1e-6;
    lambda = 1000;
    U = 200;
    OutIter = U*lambda;
    t1 = cputime;
    [Disc1, Disc1vec,ID1] = TA_REM_CD_f(n,q,OutIter,U,T0,T1,Reps,writedata);
    fprintf('%.8f\n',Disc1)


