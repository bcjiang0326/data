function Data_Lrg()
% 20130102 
% 将 试验次数在 100 以内的最优无重复试验数据整理
% 首先将数据全部读进来
datafile = fopen('LargeNoReps20121224.txt','r');
data = fscanf(datafile,'%f',[1,inf]);
fclose(datafile);

% 设置输出文件
outfile = fopen('disc.txt','w');
fprintf(outfile,'q s n & CD & q s n & CD & q s n & CD \\\\ \\hline \n');

N = length(data);
bgn = 1; %读数据开始位置
cs = 1;
while bgn < N
    levels = data(bgn);
    s = data(bgn+1);
    n = data(bgn+2);
    CD = data(bgn+3);
    out=[levels,s,n,CD];
    fprintf(outfile,'%d %d %d & %.8f',out);
    if mod(cs,3)==0
        fprintf(outfile,'\\\\ \n');
    else
        fprintf(outfile,' & ');
    end
    
    bgn = bgn+4+n;
    cs = cs + 1;
end

fclose(outfile);

end