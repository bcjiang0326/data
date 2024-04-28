Reps = 30;
dims = [2,3];
outfile = fopen('MA LHD 20151015/comp_ma36_20151016.txt','w');
for n = 3:9
    %Mm0 = zeros(Reps,1);
    %MProj0 = zeros(Reps,1);
    CD0 = zeros(Reps,1);
    ptn0 = zeros(Reps,length(dims));
    Es0 = zeros(Reps,1);
    Mmp0 = zeros(Reps,2);
    for rep = 1:Reps
        Lname = strcat('MA LHD 20151015/L0_',int2str(n),'_',int2str(rep));
        L = importdata(Lname);
        N = size(L,1);
        %Mm0(rep) = MaxminValue(L,2);
        %MProj0(rep) = MaxProjection(L);
        CD0(rep) = CD2_value(L,N*ones(n,1));
        ptn0(rep,:) = CD2_pattern(L,N*ones(n,1),dims);
        cor0 = corr(L);Es0(rep) = (sum(sum(cor0))-n)/(n*(n-1));
        Mmp0(rep,:) = Maxminproj(L,dims);
    end
    
    %Mm2 = zeros(Reps,1);
    %MProj2 = zeros(Reps,1);
    CD2 = zeros(Reps,1);
    ptn2 = zeros(Reps,length(dims));
    Es2 = zeros(Reps,1);
    Mmp2 = zeros(Reps,2);
    for rep = 1:Reps
        Lname = strcat('MA LHD 20151015/L2_',int2str(n),'_',int2str(rep));
        L = importdata(Lname);
        N = size(L,1);
        %Mm2(rep) = MaxminValue(L,2);
        %MProj2(rep) = MaxProjection(L);
        CD2(rep) = CD2_value(L,N*ones(n,1));
        ptn2(rep,:) = CD2_pattern(L,N*ones(n,1),dims);  
        cor2 = corr(L);Es2(rep) = (sum(sum(cor2))-n)/(n*(n-1));
        Mmp2(rep,:) = Maxminproj(L,dims) ;
    end
    y = [CD0,ptn0,Mmp0,CD2,ptn2,Mmp2];
    %fprintf(outfile,'%.4e %.4e %.4e %.4f %.4e %.4e %.4e %.4f\n',mean(y));
    %fprintf('%.6f %.6f %.6f %.4f %.6f %.6f %.6f %.4f\n',mean(y));
    fprintf('%.4e ',mean(y)); fprintf('\n');
end
fclose(outfile);