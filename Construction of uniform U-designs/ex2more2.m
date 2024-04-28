% 本程序去掉 TA 算法中的最小和最大MD对应的设计，换入两个新的设计，
% 使得其最小MD增大，平均 MD 保持不变

ex2more;

epsilon = 1e-10; T0 = 1e-2;  T1 = 1e-6;
InIter = 1e4; OutIter = 100*InIter; Reps = 30;

N = 36; n = 7;
%{
q_lh = N*ones(n,1);
Disc2 = zeros(Reps,1);
L2 = cell(Reps,1);
for rep = 1:Reps
    [Disc2(rep),~, L2{rep}] = TA_MD(N,q_lh,OutIter*1.4,InIter*1.4,T0,T1,1,0);
end
%}


[~,id1] = sort(Disc1);
aveDisc1 = mean(Disc1);
[Disc2,a] = sort(Disc2);
L2 = L2(a);
%
diff  = inf;
i0 = 0; j0 = 0;
for i = 1:Reps-1
    if Disc2(i)<1.0260e-1
        continue
    end
    for j = i+1:Reps
        diff_ = abs(mean([Disc2(i);Disc1(id1(2:end-1));Disc2(j)])-aveDisc1);
        if diff_ < diff-epsilon
            diff = diff_;
            i0 = i;
            j0 = j;
        end
        if diff < 1e-6
            break;
        end
    end
end
%}

%{
i0 = 0;
for i = 1:Reps-1
    if Disc2(i) > 1.0260e-1
        i0 = i;
        break;
    end
end

for j = i0+1:Reps
    diff_ = abs(mean([Disc2(i0);Disc1(id1(2:end-1));Disc2(j)])-aveDisc1);
    if diff_ < diff-epsilon
        diff = diff_;
        j0 = j;
    end
end
%}


Lname =strcat('MA LHD 20151024 MD oa36 mwa ta/L_ta_',int2str(n),'_',int2str(id1(1)));
out= fopen(Lname,'w');
for i = 1:N
    fprintf(out,'%d ',L2{i0}(i,:));
    fprintf(out,'\n');
end
fclose(out);

Lname =strcat('MA LHD 20151024 MD oa36 mwa ta/L_ta_',int2str(n),'_',int2str(id1(end)));
out= fopen(Lname,'w');
for i = 1:N
    fprintf(out,'%d ',L2{j0}(i,:));
    fprintf(out,'\n');
end
fclose(out);

ex2more;

Lname =strcat('MA LHD 20151024 MD oa36 mwa ta/L_ta_',int2str(n),'_',int2str(18));
data = importdata(Lname);
fprintf('%.8e\n',MD2_value(data,N*ones(n,1)));








