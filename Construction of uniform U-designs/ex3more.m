% 分别基于MWA准则前Reps最优和GMA准则前Reps最优来构造uniform U-design，比较二者之前的关系

epsilon = 1e-10; T0 = 1e-2;  T1 = 1e-6; 
InIter = 1e4; OutIter = 100*InIter; Reps = 100;
%InIter = 1e3; OutIter = 10*InIter; Reps = 10;

oa = importdata('OA from web/MA.36.3.12.2.11.txt');
[N,ncols] = size(oa);
q = zeros(ncols,1);
for i = 1:ncols
    q(i) = length(unique(oa(:,i)));
end
[q,id] = sort(q);
oa = oa(:,id);
n1 = sum(q==2);
n2 = sum(q==3);

nn = [2 3; 3 3; 2 4; 4 3; 3 4; 4 4];
%fid = fopen('ex3.txt','w');
for lp = 1:size(nn,1)
    n10 = nn(lp,1);
    n20 = nn(lp,2);
    n = n10+n20;
    q_oa = [2*ones(n10,1); 3*ones(n20,1)];
    filename = 'MA LHD 20151102 MD oa36 comp order/';
    filename = strcat(filename,'results',int2str(n10),int2str(n20));
    
    results = importdata(filename);
    cols = results(:,1:n);
    ia = results(:,n+3);
    [~,ia] = sort(ia);
    iwb = results(:,n+4);
    [~,iwb] = sort(iwb);
    clear('results');
    
    %{
    Disc_vec = zeros(size(cols,1),1);    
    for i = 1:size(cols,1)
        if ia(i) > Reps && iwb(i) > Reps
            continue;
        else
            [Disc_vec(i), ~, ~] = TA_MD_LHD(oa(:,cols(i,:)),q_oa,OutIter,InIter,T0,T1,1,0,1);
        end
    end
    out = fopen(strcat(filename,'Disc'),'w');
    fprintf(out,'%.8e ',Disc_vec);
    fclose(out);
    %}
    id1 = iwb<=Reps;
    id2 = ia<=Reps;
    fprintf('%d\n',sum(id1~=id2));
    %y = [n10,n20,mean(Disc_vec(id1)),std(Disc_vec(id1)),min(Disc_vec(id1)),...
    %    mean(Disc_vec(id2)),std(Disc_vec(id2)),min(Disc_vec(id2))];
    %fprintf(fid,'%d %d %.4e %.4e %.4e %.4e %.4e %.4e\n',y);
    %fprintf('%d %d  %.4e %.4e %.4e  %.4e %.4e %.4e\n',y);    
end
%fclose(fid);