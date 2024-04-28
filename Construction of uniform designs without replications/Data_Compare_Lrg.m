function Data_Compare_Lrg( )
% 20130106
% 比较若干文件的数据值，优选 CD 小的，其次选 WD 小的
datafile1 = fopen('LargeNoReps20121224.txt','r');
datafile2 = fopen('LargeNoReps20130107.txt','r');
data1 = fscanf(datafile1,'%f',[1,inf]);
data2 = fscanf(datafile2,'%f',[1,inf]);
fclose(datafile1); fclose(datafile2);

% 设置输出文件
out = fopen('LargeNoReps20130109.txt','w');

epsilon = 1e-12;
N = length(data1);
bgn = 1; %读数据开始位置
%cs = 1; %case number
while bgn < N
    levels = data1(bgn);
    s = data1(bgn+1);
    n = data1(bgn+2);
    q = levels*ones(s,1);
    CD1 = data1(bgn+3);
    CD2 = data2(bgn+3);
    if CD1 < CD2
        CD = CD1;
        y = data1(bgn+4:bgn+3+n)';
        cout_d([levels,s,n,2]);
    elseif CD1 > CD2
        CD = CD2;
        y = data2(bgn+4:bgn+3+n)';
        cout_d([levels,s,n,2]);
    else
        CD = CD1;
        y1 = data1(bgn+4:bgn+3+n)';
        y2 = data2(bgn+4:bgn+3+n)';
        D1 = Id2Design(q,y1);
        D2 = Id2Design(q,y2);
        WD1 = WD2_value(D1,q);
        WD2 = WD2_value(D2,q);
        if WD1 < WD2-epsilon
            y = y1;
            cout_d([levels,s,n,1]);
        else
            y = y2;
            cout_d([levels,s,n,2]);
        end
    end
    
    fprintf(out,'%d %d %d %.8f ',[levels,s,n,CD]);
    fprintf(out,'%d ',y);
    fprintf(out,'\n');
    bgn = bgn+4+n;
end

fclose(out);

