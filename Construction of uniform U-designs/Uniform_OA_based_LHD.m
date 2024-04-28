function Uniform_OA_based_LHD( )
data = importdata('uniform OA/CD_results.txt');
%data = importdata('uniform OA/CD_results3.txt');
K = size(data,1);

outfile = fopen('uniform LHD/CD_results3.txt','w');
for lp = 1:K
    Disc = inf;
    N = data(lp,1); s = data(lp,2); n = data(lp,3);
    q = s*ones(n,1);
    D0 = importdata(['uniform OA/N',int2str(N),'s',int2str(s),'n',int2str(n),'.txt']);
    for k = 1:30
        InIter = 200;
        OutIter = 1000*InIter;
        Reps = 30;
        [Disc0, Disc_vec0, L0] = TA_CD_LHD(D0,q,OutIter,InIter,1e-2,1e-6,Reps,0);
        %fprintf('%.8f\n',Disc0);
        if Disc0<Disc-1e-10
            Disc = Disc0;
            Disc_vec = Disc_vec0;
            L = L0;
        end
    end
    L = sortrows(L); 
    
    fprintf(outfile,'%d %d %.6e %.6e\n',[N,n,mean(Disc_vec),Disc]);
    fprintf('%d %d %.6e %.6e\n',[N,n,mean(Disc_vec),Disc]);
    
    str = ['uniform LHD/',int2str(N),'_',int2str(n),'_s',int2str(s),'.txt'];
    fid = fopen(str,'w');
    for i  = 1:N
        fprintf(fid,'%d ',L(i,:));
        fprintf(fid,'\n');
    end
    fclose(fid);
end
fclose(outfile);