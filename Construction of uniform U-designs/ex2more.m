%¼ÆËã30¸öUDµÄmean std£¬min 

OA = importdata('OA from web/MA.36.3.12.2.11.txt');
[N,ncols] = size(OA);
q = zeros(ncols,1);
for i = 1:ncols
    q(i) = length(unique(OA(:,i)));
end
[q,id] = sort(q);
OA = OA(:,id);

Reps = 30;
array = [2,3; 
    3,3; 
    2,4; 
    4,3; 
    3,4; 
    4,4];
%array = [3,4];
DiscFunc = @MD2_value;
for i = 1:size(array,1)
    n10 = array(i,1);
    n20 = array(i,2);
    n = n10+n20;
    q_lh = N*ones(n,1);
    Disc0 = zeros(Reps,1);
    for rep = 1: Reps
        Lname = strcat('MA LHD 20151024 MD oa36 mwa ta/L_mwa',int2str(n10),int2str(n20),'_',int2str(rep));
        L = importdata(Lname);
        Disc0(rep) = DiscFunc(L,q_lh);
    end
    fprintf('%d %d: %.4e %.4e %.4e\n',n10,n20,mean(Disc0),std(Disc0),min(Disc0));
end

narray = unique(sum(array,2));
for i = 1:size(narray,1)
    n = narray(i);
    q_lh = N*ones(n,1);
    Disc1 = zeros(Reps,1);
    for rep = 1: Reps
        Lname = strcat('MA LHD 20151024 MD oa36 mwa ta/L_ta_',int2str(n),'_',int2str(rep));
        L = importdata(Lname);
        %if rep==12
        %    a = 4;
        %end
        Disc1(rep) = DiscFunc(L,q_lh);
    end
    fprintf('%d: %.4e %.4e %.4e\n',n,mean(Disc1),std(Disc1),min(Disc1));
end