function uniform_OA( )
allname = {'OA from web/oa.16.5.4.2.a.txt',...
    'OA from web/oa.16.8.2.3.txt',...
    'OA from web/oa.18.7.3.2.txt',...
    'OA from web/oa.24.12.2.3.txt',...
    'OA from web/oa.25.6.5.2.txt',...
    'OA from web/oa.32.9.4.2.a.txt',...
    'OA from web/oa.36.13.3.2.txt',...
    'OA from web/oa.50.11.5.2.txt'};

fid0 = fopen('uniform OA/CD_results.txt','w');
for k = 1:length(allname)
    OA = importdata(allname{k});
    [N,n1] = size(OA); s = length(unique(OA(:,1))); 
    for n = 3:n1
        name = ['uniform OA/N',int2str(N),'s',int2str(s),'n',int2str(n),'.txt'];
        fid = fopen(name,'w');
        
        D = MA_from_OA(OA,n);
    
        InIter = 200;
        OutIter = 1000*InIter;
        q = s*ones(n,1);
        [CD0, ~, D0] = TA_uniform_FF_design(D,q,OutIter,InIter,1e-2,1e-5,50,0);
        D0 = sortrows(D0);
        
        fprintf('%d %d %d %.6f\n',[N,s,n,CD0]);
        fprintf(fid0,'%d %d %d %.6f\n',[N,s,n,CD0]);    
        for i = 1:N
            fprintf(fid,'%d ',D0(i,:));
            fprintf(fid,'\n');
        end
        fclose(fid);
    end
end
fclose(fid0);


    