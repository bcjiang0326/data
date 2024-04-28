function Data_Smll()
% 20130101 
% 处理网站上所有有重复的设计
% 首先将数据全部读进来
datafile = fopen('SmallNoReps20130106.txt','r');
data = fscanf(datafile,'%f',[1,inf]);
fclose(datafile);

% 设置输出文件
outf1 = fopen('disc.txt','w');
outf2 = fopen('design.txt','w');
fprintf(outf1,'q s n & CD & WD & q s n & CD & WD \\\\ \\hline \n');
fprintf(outf2,'q s n & f(D) \\\\ \\hline \n');
%div_ave = zeros(2,1);
epsilon = 1e-12;
N = length(data);
bgn = 1; %读数据开始位置
cs = 0; %case number
while bgn < N
    cs = cs+1;
    
    levels = data(bgn);
    s = data(bgn+1);
    n = data(bgn+2);
    q = levels*ones(s,1);
    y = data(bgn+4:bgn+3+n)';
    if levels == 3
        fname = strcat('CD2/Level 3/',int2str(s),'_',int2str(n),'.txt');
    elseif levels == 4
        fname = strcat('CD2/Level 4/',int2str(s),'_',int2str(n),'.txt');
    elseif levels == 5
        fname = strcat('CD2/Level 5/',int2str(s),'_',int2str(n),'.txt');
    elseif levels == 6
        fname = strcat('CD2/Level 6/',int2str(s),'_',int2str(n),'.txt');
    end
    D = Id2Design(q,y);
    Dw = importdata(fname);Dw = Dw-1;
    CD = CD2_value(D,q);
    CDw = CD2_value(Dw,q);
    div_CD = 100*(CD-CDw)/CDw;
    WD = WD2_value(D,q);
    WDw = WD2_value(Dw,q);
    div_WD = 100*(WD-WDw)/WDw;
    if abs(div_CD) < epsilon
        div_CD = 0;
    end
    if abs(div_WD) < epsilon
        div_WD = 0;
    end
%    div_ave = div_ave + [div_CD;div_WD];
    out=[levels,s,n,CD,div_CD,WD,div_WD];
    fprintf(outf1,'%d %d %d & %.6f(%.2f) & %.6f(%.2f)',out);
    if mod(cs,2)==0
        fprintf(outf1,'\\\\ \n');
    else
        fprintf(outf1,' & ');
    end
    
    fprintf(outf2,'%d %d %d &',[levels,s,n]);
    fprintf(outf2,'%d ', y);
    fprintf(outf2,'\\\\ \n');
    
    bgn = bgn+4+n;
end

%div_ave = div_ave/cs;
%cout_8f(div_ave); cout_d(cs);
fclose(outf1);fclose(outf2);

end