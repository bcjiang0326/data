
name = {
    'OA from web/oa.16.15.2.2.0.txt',...
    'OA from web/oa.24.12.2.3.txt'};
%epsilon = 1e-10; InIter = 10^4;  OutIter = 100*InIter;
%T0 = 1e-2;  T1 = 1e-6; Reps = 100;
fprintf('CD:\n');
outfile = fopen('MA LHD 20150126/result_CD_20150525.txt','w');
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
        groups = zeros(size(D)); 
        for j = 1:n
            [~,groups(:,j)] = sort(D(:,j));
        end
        t0 = cputime;
        [Disc0, Disc_vec0] = TA_CD_LHD_1(groups,D,q_oa,OutIter,InIter,T0,T1,Reps,0);
        t0 = cputime-t0;


        % 2.计算 uniform OA-based LHD
        t1 = 0;
        Disc_vec1 = zeros(Reps,1);
        Disc1 = inf;
        for rep = 1:Reps
            id = sort(randperm(ncols,n));
            D = OA(:,id);
            groups = zeros(size(D)); 
            for j = 1:n
                [~,groups(:,j)] = sort(D(:,j));
            end
            tt = cputime;
            Disc_vec1(rep) = TA_CD_LHD_1(groups,D,q_oa,OutIter,InIter,T0,T1,1,0);
            if Disc_vec1(rep)<Disc1-epsilon
                Disc1 = Disc_vec1(rep);
                %L1 = L;
            end
            tt = cputime-tt;
            t1 = t1+tt;
        end
        
        % 3.直接采用 TA 算法 计算 uniform LHD
        t2 = cputime;
        [Disc2, Disc_vec2, L2,K] = TA_CD(N,q_lh,OutIter,InIter,T0,T1,Reps,0);
        t2 = cputime-t2;
        
        y = [N,n,mean(Disc_vec0),t0,mean(Disc_vec1),t1,mean(Disc_vec2),t2];
        fprintf(outfile,'%d  &  %d  &  %.6f  &  %.2f  &  %.6f  &  %.2f  &  %.6f  &  %.2f  \\\\ \n',y);
        fprintf('%d  &  %d  &  %.6f  &  %.2f  &  %.6f  &  %.2f  &  %.6f  &  %.2f  \\\\ \n',y);
    end
end
fclose(outfile);