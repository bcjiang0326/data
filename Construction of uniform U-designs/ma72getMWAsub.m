
OA = importdata('OA from web/MA.72.4.1.3.24.2.20.txt');
[N,ncols] = size(OA);
q = zeros(ncols,1);
for i = 1:ncols
    q(i) = length(unique(OA(:,i)));
end
[q,id] = sort(q);
OA = OA(:,id);

lv = [2,3,4];



M = 10;
eps = 1e-12;

for k = 2:6
    n1 = k;
    n2 = k;
    n3 = 1;
    aveDisc0 = inf;
    v0 = [];
    Disc_vec = zeros(M,1);
    for i = 1:M
        cols1 = randperm(20,n1);
        cols2 = randperm(24,n2)+20;
        v = sort([cols1,cols2,45]);
        D = OA(:,v);
        q0 = [2*ones(n1,1);3*ones(n2,1);4];
        [~,aveDisc] = Weighted_wordtype_pattern(D,q0,'MD');
        Disc_vec(i) = aveDisc;
        if aveDisc < aveDisc0-eps
            aveDisc0 = aveDisc;
            v0 = v;
        end
    end
    y = [n1,n2,n3,aveDisc0,mean(Disc_vec)];
    fprintf('%d %d %d %.6f %.6f:   ',y)
    fprintf('%d ',v)
    fprintf('\n')
end