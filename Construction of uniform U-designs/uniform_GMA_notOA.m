function uniform_GMA_notOA( )
% 20131022
% ËÑË÷ uniform GMA design
% ËÑË÷ uniform GMA design based LHD

fid0 = fopen('uniform GMA/CD_results2.txt','a');
fid1 = fopen('uniform LHD1/CD_results2.txt','a');

Nzone = 26; nzone = [8,9,10]; s = 2;
for N = Nzone;
    for n = nzone
        Iter = 1e5;
        Reps = 100;
        [~, ~, D] = GMA(N,n,s,Iter,Reps);
        
        InIter = 200;
        OutIter = 1000*InIter;
        q = s*ones(n,1);
        if s~=2
            [CD0, ~, D0] = TA_uniform_FF_design(D,q,OutIter,InIter,1e-2,1e-5,1000,0);
            D0 = sortrows(D0);
        else
            D0 = sortrows(D);
            CD0 = CD2_value(D0,q);
        end
        
        name = ['uniform GMA/N',int2str(N),'s',int2str(s),'n',int2str(n),'.txt'];
        outUGMA = fopen(name,'w');
        fprintf('UGMA: %d  %d  %d  %.6f\n',[N,s,n,CD0]);
        fprintf(fid0,'%d %d %d %.6f\n',[N,s,n,CD0]);    
        for i = 1:N
            fprintf(outUGMA,'%d ',D0(i,:));
            fprintf(outUGMA,'\n');
        end
        fclose(outUGMA);
        
        Disc = inf;
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
        
        fprintf('LHD: %d  %d  %.6e  %.6e\n',[N,n,mean(Disc_vec),Disc]);
        fprintf(fid1,'%d %d %.6e %.6e\n',[N,n,mean(Disc_vec),Disc]);
        
        str = ['uniform LHD1/',int2str(N),'_',int2str(n),'_s',int2str(s),'.txt'];
        outLHD = fopen(str,'w');
        for i  = 1:N
            fprintf(outLHD,'%d ',L(i,:));
            fprintf(outLHD,'\n');
        end
        fclose(outLHD);
    end
end

fclose(fid0); fclose(fid1);
        
        
        
        