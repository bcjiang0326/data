InIter = 500;
OutIter = 200*InIter;
Reps = 1;
dims = 2;
epsilon = 1e-12;

s = 4;
for n = 4:15
    for N = 8:s:24
        fname = strcat('../CD2/Level 4/',int2str(n),'_',int2str(N),'.txt');
        D = importdata(fname); D = D-1;
        q = s*ones(n,1);
        ctn = CD2_pattern(D,q,2:n);
        ptn = CD2_pattern(D,q,dims);
        [ptn0, ptn_vec, D0] = TA_proj_CD(dims, N,q,OutIter,InIter,1e-2,1e-6,Reps,0);
        CD = CD2_value(D,q);
        CD0 = CD2_value(D0,q);
        v = find(abs(ptn0-ptn)>epsilon,1);
        if ptn0(v)<ptn(v)
             fprintf('n = %d, N = %d, delta = %.4e, deltaCD = %.4e\n',[n,N,ptn0-ptn,CD0-CD]);
             str = strcat('results/Level',int2str(s),'_N',int2str(N),'n',int2str(n),'.txt');
             fid = fopen(str,'w');
             for i = 1:N
                 fprintf(fid,'%d ',D0(i,:));
                 fprintf(fid,'\n');
             end
             fclose(fid);
        end
    end
end

s = 5;
for n = 4:12
    for N = 10:s:25
        fname = strcat('../CD2/Level 5/',int2str(n),'_',int2str(N),'.txt');
        D = importdata(fname); D = D-1;
        q = s*ones(n,1);
        ctn = CD2_pattern(D,q,2:n);
        ptn = CD2_pattern(D,q,dims);
        [ptn0, ptn_vec, D0] = TA_proj_CD(dims, N,q,OutIter,InIter,1e-2,1e-6,Reps,0);
        CD = CD2_value(D,q);
        CD0 = CD2_value(D0,q);
        v = find(abs(ptn0-ptn)>epsilon,1);
        if ptn0(v)<ptn(v)
             fprintf('n = %d, N = %d, delta = %.4e, deltaCD = %.4e\n',[n,N,ptn0-ptn,CD0-CD]);
             str = strcat('results/Level',int2str(s),'_N',int2str(N),'n',int2str(n),'.txt');
             fid = fopen(str,'w');
             for i = 1:N
                 fprintf(fid,'%d ',D0(i,:));
                 fprintf(fid,'\n');
             end
             fclose(fid);
        end
    end
end

s = 6;
for n = 4:9
    for N = 12:s:30
        fname = strcat('../CD2/Level 6/',int2str(n),'_',int2str(N),'.txt');
        D = importdata(fname); D = D-1;
        q = s*ones(n,1);
        ctn = CD2_pattern(D,q,2:n);
        ptn = CD2_pattern(D,q,dims);
        [ptn0, ptn_vec, D0] = TA_proj_CD(dims, N,q,OutIter,InIter,1e-2,1e-6,Reps,0);
        CD = CD2_value(D,q);
        CD0 = CD2_value(D0,q);
        if ptn0<ptn-1e-10
             fprintf('n = %d, N = %d, delta = %.4e, deltaCD = %.4e\n',[n,N,ptn0-ptn,CD0-CD]);
             str = strcat('results/Level',int2str(s),'_N',int2str(N),'n',int2str(n),'.txt');
             fid = fopen(str,'w');
             for i = 1:N
                 fprintf(fid,'%d ',D0(i,:));
                 fprintf(fid,'\n');
             end
             fclose(fid);
        end
    end
end


s = 3;
for n = 4:15
    for N = 9:s:24
        fname = strcat('../CD2/Level 3/',int2str(n),'_',int2str(N),'.txt');
        D = importdata(fname); D = D-1;
        q = s*ones(n,1);
        ctn = CD2_pattern(D,q,2:n);
        ptn = CD2_pattern(D,q,dims);
        [ptn0, ptn_vec, D0] = TA_proj_CD(dims, N,q,OutIter,InIter,1e-2,1e-6,Reps,0);
        CD = CD2_value(D,q);
        CD0 = CD2_value(D0,q);
        if ptn0<ptn-1e-10
             fprintf('n = %d, N = %d, delta = %.4e, deltaCD = %.4e\n',[n,N,ptn0-ptn,CD0-CD]);
             str = strcat('results/Level',int2str(s),'_N',int2str(N),'n',int2str(n),'.txt');
             fid = fopen(str,'w');
             for i = 1:N
                 fprintf(fid,'%d ',D0(i,:));
                 fprintf(fid,'\n');
             end
             fclose(fid);
        end
    end
end

