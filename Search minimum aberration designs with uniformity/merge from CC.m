
OAname = {'OA from web/oa.25.6.5.2.txt';
    'OA from web/oa.32.9.4.2.a.txt';
    'OA from web/oa.36.13.3.2.txt';
    'OA from web/oa.36.3.6.2.txt';
    'OA from web/oa.50.11.5.2.txt'};
n_ed = [12,16,15,12,12];
n_bg = [12,13,15,6,12];

for k = [2,4]
    oa = importdata(OAname{k});
    [N,n_oa] = size(oa);
    s = length(unique(oa(:,1)));
    for n = n_bg(k)+1:n_ed(k)
        outfile1 = strcat('MAD.',int2str(N),'.',int2str(s),'/MAD.',int2str(N),'.',int2str(s),'.',int2str(n),'.txt');
        outfile2 = strcat('From CC/MAD.',int2str(N),'.',int2str(s),'/MAD.',int2str(N),'.',int2str(s),'.',int2str(n),'.txt');
        
        D1 = importdata(outfile1);
        [aveCD1,A1] = aveCD_LevelPerm(D1);
        
        D2 = importdata(outfile2);
        [aveCD2,A2] = aveCD_LevelPerm(D2);
        
        %
        if aveCD2 < aveCD1-1e-10
            fid = fopen(outfile1,'w');
            for i = 1:N
                fprintf(fid,'%d ',D2(i,:));
                fprintf(fid,'\n');
            end
            fclose(fid);
        end
        

        fprintf('%d %d %d: %.6f  (%.2f,%.2f,%.2f)  %.6f  (%.2f,%.2f,%.2f)\n', N,s,n,aveCD1,A1(2:4),aveCD2,A2(2:4));
    end
end
