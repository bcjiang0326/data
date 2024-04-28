Reps = 30;
dims = [2,3];
N = 64;
%outfile = fopen('MALHD 20151020 oa64/comp_oa64_20151021.txt','w');
%zone = [3:20];
zone = [4];
for n = zone
    q_lh = N*ones(n,1);
    MD0 = zeros(Reps,1);
    ptn0 = zeros(Reps,length(dims));
    %Mmp0 = zeros(Reps,length(dims));
    %Es0 = zeros(Reps,1);
    for rep = 1:Reps
        Lname = strcat('MALHD 20151020 oa64/L0_',int2str(N),'_',int2str(n),'_',int2str(rep));
        L = importdata(Lname);
        N = size(L,1);
        MD0(rep) = MD2_value(L,q_lh);
        ptn0(rep,:) = MD2_pattern(L,q_lh,dims);
        %Mmp0(rep,:) = Maxminproj(L,dims);
        %cor0 = corr(L);Es0(rep) = (sum(sum(cor0))-n)/(n*(n-1));
    end
    [~,id] = sort(MD0);
    MD0(id(end),:) = [];
    ptn0(id(end),:) = [];
    
    MD2 = zeros(Reps,1);
    ptn2 = zeros(Reps,length(dims));
    %Mmp2 = zeros(Reps,length(dims));
    %Es2 = zeros(Reps,1);
    for rep = 1:Reps
        Lname = strcat('MALHD 20151020 oa64/L2_',int2str(N),'_',int2str(n),'_',int2str(rep));
        L = importdata(Lname);
        N = size(L,1);
        MD2(rep) = MD2_value(L,q_lh);
        ptn2(rep,:) = MD2_pattern(L,q_lh,dims);  
        %Mmp2(rep,:) = Maxminproj(L,dims) ;
        %cor2 = corr(L);Es2(rep) = (sum(sum(cor2))-n)/(n*(n-1));
    end
    y = [mean(MD0),mean(ptn0),mean(MD2),mean(ptn2)];
    %fprintf(outfile,'%.4e %.4e %.4e %.4f %.4e %.4e %.4e %.4f\n',mean(y));
    %fprintf('%.6f %.6f %.6f %.4f %.6f %.6f %.6f %.4f\n',mean(y));
    fprintf('%d & ',n); fprintf('%.3e & ',y); fprintf('\\\\ \n');
end
%fclose(outfile);