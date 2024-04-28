function Uniform_OA_based_LHD_27run( )
N = 27; 
X = combins(3,3); % втсиап
C = [1 1 1; 1 2 0; 1 1 2; 1 0 1; 0 1 2; 1 2 2; 1 1 0]';
B = { [0,1], [0,1,2], [1,1,0,2], [1,1,0,2,1], [1,1,0,2,1,2], [1,1,0,2,1,2,2] };


y = [3.202036e-02, 5.094166e-02, 7.597447e-02];

outfile = fopen('uniform LHD/CD_results27.txt','w');
for n = 8:10
    Disc = inf;
    q = 3*ones(n,1);
    D0 = mod( [X,X*C(:,1:n-3)+ones(N,1)*B{n-4}], 3);
    for k = 1:200
        InIter = 200;
        OutIter = 1000*InIter;
        Reps = 30;
        [Disc0, Disc_vec0, L0] = TA_CD_LHD(D0,q,OutIter,InIter,1e-2,1e-6,Reps,0);
        if Disc0<Disc-1e-10
            Disc = Disc0;
            Disc_vec = Disc_vec0;
            L = L0;
        end
    end
    L = sortrows(L); 
 
    if Disc < y(n-7)-1e-10
        fprintf(outfile,'%d %d %.6e %.6e\n',[N,n,mean(Disc_vec),Disc]);
        fprintf('%d %d %.6e %.6e\n',[N,n,mean(Disc_vec),Disc]);
    
        str = ['uniform LHD/',int2str(N),'_',int2str(n),'.txt'];
        fid = fopen(str,'w');
        for i  = 1:N
            fprintf(fid,'%d ',L(i,:));
            fprintf(fid,'\n');
        end
        fclose(fid);
    end
end
fclose(outfile);