function checkReps_Level6( )
outfile1 = fopen('CD2/Level 6/WithReps.txt','w');
outfile2 = fopen('CD2/Level 6/Reps.txt','w');
fprintf(outfile2,'s n Point \t\t ID Nreps\n');
outfile3 = fopen('CD2/Level 6/info.txt','w');
fprintf(outfile3,'s n ID Nreps\n');

levels = 6;
for s = 2:5
    for n = 12:levels:48
        if levels^s > n
            fname = strcat('CD2/Level 6/',int2str(s),'_',int2str(n),'.txt');
            D = importdata(fname); D = D-1;
            q = levels*ones(s,1);
            [Points, NumReps, ID] = checkRep(D,q);
            if ~isempty(Points)
                fprintf(outfile1,'%d %d\n',s,n);
                for i = 1:length(ID)
                    fprintf(outfile2,'%d %d ',s,n);
                    for j = 1:s
                        fprintf(outfile2,'%d ',Points(i,j));
                    end
                    fprintf(outfile2,'%d %d\n',ID(i),NumReps(i));
                end
                for i = 1:length(ID)
                    fprintf(outfile3,'%d %d %d %d\n',s,n,ID(i),NumReps(i));
                end
            end
        end
    end
end

for s = 6:9
    for n = 12:levels:60
        if levels^s > n
            fname = strcat('CD2/Level 6/',int2str(s),'_',int2str(n),'.txt');
            D = importdata(fname); D = D-1;
            q = levels*ones(s,1);
            [Points, NumReps, ID] = checkRep(D,q);
            if ~isempty(Points)
                fprintf(outfile1,'%d %d\n',s,n);
                for i = 1:length(ID)
                    fprintf(outfile2,'%d %d ',s,n);
                    for j = 1:s
                        fprintf(outfile2,'%d ',Points(i,j));
                    end
                    fprintf(outfile2,'%d %d\n',ID(i),NumReps(i));
                end
                for i = 1:length(ID)
                    fprintf(outfile3,'%d %d %d %d\n',s,n,ID(i),NumReps(i));
                end    
            end
        end
    end
end

end

