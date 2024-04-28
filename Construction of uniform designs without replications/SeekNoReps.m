function SeekNoReps( )
% 搜索无重复设计，更新网站上面所有有重复的设计
% 并将搜索结果输入至文件
data = [3 3 18; 3 3 21; 3 3 24; ...
    3 4 15; 3 4 24; 3 4 33; 3 4 36; 3 4 42; 3 4 45; 3 4 48; 3 4 51;...
    3 5 15; 3 5 24; 3 5 33; 3 5 36; 3 5 45; 3 5 48; 3 5 51;...
    3 6 15; 3 6 24; 3 6 33; 3 6 42; 3 6 51;...
    3 7 24; 3 7 33; 3 7 42; 3 7 51;...
    3 8 33; 3 8 42; 3 8 51;...
    3 9 51;...
    4 2 12;...
    4 3 28; 4 3 36; 4 3 40; 4 3 44; 4 3 48; 4 3 52;...
    5 2 10; 5 2 15; 5 2 20;...
    5 3 10; 5 3 35;...
    5 4 35;...
    5 5 35;...
    6 2 24; 6 2 30];

N = size(data,1);
Reps = 100;
Iter = 100000;
InIter = 200;
T0 = 0.01;
T1 = 1e-6;
epsilon = 1e-11;

outfile = fopen('SmallNoReps.txt','w');
for i = 1:N
    levels = data(i,1);
    s = data(i,2); n = data(i,3);
    cout_d([levels,s,n]);
    q = levels*ones(s,1);
    Disc0 = inf;
    y0 = 0;
    for j = 1:Reps
        [Disc, y] = TA_REM_CD(n,q,Iter,InIter,T0,T1);
        if Disc < Disc0-epsilon
            Disc0 = Disc;
            y0 = y;
        end
    end
    fprintf(outfile,'%d %d %d %.8f ',levels,s,n,Disc0);
    fprintf(outfile,'%d ',y0);
    fprintf(outfile,'\n');
end
fclose(outfile);


outfile2 = fopen('LargeNoReps.txt','w');

levels = 3; 
for s = 5:15
    for n = 54:levels:99
        cout_d([levels,s,n]);
        q = levels*ones(s,1);
        Disc0 = inf;
        y0 = 0;
        for j = 1:Reps
            [Disc, y] = TA_REM_CD(n,q,Iter,InIter,T0,T1);
            if Disc < Disc0-epsilon
                Disc0 = Disc;
                y0 = y;
            end
        end
        fprintf(outfile2,'%d %d %d %.8f ',levels,s,n,Disc0);
        fprintf(outfile2,'%d ',y0);
        fprintf(outfile2,'\n');
    end
end

levels = 4;
for s = 5:15
    for n = 56:levels:100
        cout_d([levels,s,n]);
        q = levels*ones(s,1);
        Disc0 = inf;
        y0 = 0;
        for j = 1:Reps
            [Disc, y] = TA_REM_CD(n,q,Iter,InIter,T0,T1);
            if Disc < Disc0-epsilon
                Disc0 = Disc;
                y0 = y;
            end
        end
        fprintf(outfile2,'%d %d %d %.8f ',levels,s,n,Disc0);
        fprintf(outfile2,'%d ',y0);
        fprintf(outfile2,'\n');
    end
end

levels = 5;
for s = 5:12
    for n = 60:levels:100
        cout_d([levels,s,n]);
        q = levels*ones(s,1);
        Disc0 = inf;
        y0 = 0;
        for j = 1:Reps
            [Disc, y] = TA_REM_CD(n,q,Iter,InIter,T0,T1);
            if Disc < Disc0-epsilon
                Disc0 = Disc;
                y0 = y;
            end
        end
        fprintf(outfile2,'%d %d %d %.8f ',levels,s,n,Disc0);
        fprintf(outfile2,'%d ',y0);
        fprintf(outfile2,'\n');
    end
end

fclose(outfile2);

end
