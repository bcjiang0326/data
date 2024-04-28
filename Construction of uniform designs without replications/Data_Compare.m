function Data_Compare( )
% 20130106
% 比较若干文件的数据值，优选 CD 小的，其次选 WD 小的
datafile1 = fopen('SmallNoReps20121213.txt','r');
datafile2 = fopen('SmallNoReps20121216.txt','r');
datafile3 = fopen('SmallNoReps20130105.txt','r');
data1 = fscanf(datafile1,'%f',[1,inf]);
data2 = fscanf(datafile2,'%f',[1,inf]);
data3 = fscanf(datafile3,'%f',[1,inf]);
fclose(datafile1); fclose(datafile2); fclose(datafile3);

% 设置输出文件
out = fopen('SmallNoReps20130106.txt','w');

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
        CD12 = CD1;
        y12 = data1(bgn+4:bgn+3+n)';
    elseif CD1 > CD2
        CD12 = CD2;
        y12 = data2(bgn+4:bgn+3+n)';
    else
        CD12 = CD1;
        y1 = data1(bgn+4:bgn+3+n)';
        y2 = data2(bgn+4:bgn+3+n)';
        D1 = Id2Design(q,y1);
        D2 = Id2Design(q,y2);
        WD1 = WD2_value(D1,q);
        WD2 = WD2_value(D2,q);
        if WD1 < WD2-epsilon
            y12 = y1;
        else
            y12 = y2;
        end
    end
    
    CD3 = data3(bgn+3);
    y3 = data3(bgn+4:bgn+3+n)';
    if CD12 < CD3
        CD = CD12;
        y = y12;
        cout_d([levels,s,n,12]);
    elseif CD12 > CD3
        CD = CD3;
        y = y3;
        cout_d([levels,s,n,3]);
    else
        CD = CD3;
        D12 = Id2Design(q,y12);
        D3 = Id2Design(q,y3);
        WD12 = WD2_value(D12,q);
        WD3 = WD2_value(D3,q);
        if WD12 < WD3-epsilon
            y = y12;
            cout_d([levels,s,n,12]);
        else
            y = y3;
            cout_d([levels,s,n,3]);
        end
    end
    fprintf(out,'%d %d %d %.8f ',[levels,s,n,CD]);
    fprintf(out,'%d ',y);
    fprintf(out,'\n');
    bgn = bgn+4+n;
end

fclose(out);

