
name = {
    'OA from web/oa.16.15.2.2.0.txt',...
    'OA from web/oa.27.13.3.2.txt'};
%epsilon = 1e-10; InIter = 10^4;  OutIter = 100*InIter;
%T0 = 1e-2;  T1 = 1e-6; Reps = 100;
fprintf('WD:\n');
outfile = fopen('MA LHD 20150126/result_WD_20150519.txt','w');
fprintf(outfile,'InIter = %d;  OutIter = %d*InIter;  Reps = %d;\n',InIter,OutIter/InIter,Reps);
fprintf(outfile,'  &   &     MA-based LHD    &  non-MA-based LHD  &  TA  &  Existed \\\\ \n');
fprintf(outfile,'N & n & ACD & MinCD & time & AveCD & MinCD & time & ACD & MinCD & time & MinCD \\\\ \n');
for lop = 1:length(name)
    OA = importdata(name{lop});
    s = length(unique(OA(:,1)));
    [N,ncols] = size(OA);
    for n = 3:2:9
        q_oa = s*ones(n,1);
        q_lh = N*ones(n,1);
        
        % 1.计算 uniform MA-based LHD
        [D,gwp,id] = MA_from_OA(OA,n);
        t0 = cputime;
        [Disc0, Disc_vec0, L0] = TA_WD_LHD(D,q_oa,OutIter,InIter,T0,T1,Reps,0);
        t0 = cputime-t0;
        % 写入到文件中
        str = ['MA LHD 20140709/',int2str(N),'_',int2str(n),'.txt'];
        fid = fopen(str,'w');
        for i  = 1:N
            fprintf(fid,'%d ',L0(i,:));
            fprintf(fid,'\n');
        end
        fclose(fid);

        % 2.计算 uniform OA-based LHD
        t1 = cputime;
        Disc_vec1 = zeros(Reps,1);
        Disc1 = inf;
        for rep = 1:Reps
            id = sort(randperm(ncols,n));
            [Disc_vec1(rep),~, L] = TA_WD_LHD(OA(:,id),q_oa,OutIter,InIter,T0,T1,1,0);
            if Disc_vec1(rep)<Disc1-epsilon
                Disc1 = Disc_vec1(rep);
                L1 = L;
            end
        end
        t1 = cputime-t1;
        
        % 3.直接采用 TA 算法 计算 uniform LHD
        t2 = cputime;
        [Disc2, Disc_vec2, L2] = TA_WD(N,q_lh,OutIter,InIter,T0,T1,Reps,0);
        t2 = cputime-t2;
        
        y = [N,n,mean(Disc_vec0),Disc0,t0,mean(Disc_vec1),Disc1,t1,mean(Disc_vec2),Disc2,t2];
        fprintf(outfile,'%d  &  %d  &  %.6f  &  %.6f  &  %.2f  &  %.6f  &  %.6f  &  %.2f  &  %.6f  &  %.6f  &  %.2f  & \\\\ \n',y);
        fprintf('%d  &  %d  &  %.6f  &  %.6f  &  %.2f  &  %.6f  &  %.6f  &  %.2f  &  %.6f  &  %.6f  &  %.2f  & \\\\ \n',y);
    end
end
fclose(outfile);