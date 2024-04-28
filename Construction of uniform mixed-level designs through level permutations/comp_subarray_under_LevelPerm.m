%目的：对于给定 Mixed OA的所有可能子表，计算子表的 aveDisc,Cols,A,BB，
%     并将所有参数相同的子表写入一个文件.
%输出：例如 Results/MA.36.2.11.3.12.sub.3.5.txt,表示MA.36.2.11.3.12 的所有可能子表
%      MA.36.2.3.3.5 的 aveDisc,Cols,A,BB 等结果


%oa = importdata('../Web OA/MA.36.3.7.6.3.finney.txt');
oa = importdata('../Web OA/MA.36.2.11.3.12.txt');
%oa = importdata('../Web OA/MA.36.3.12.12.1.txt');
%oa = importdata('../Web OA/MA.72.3.24.24.1.txt');
q = zeros(size(oa,2),1);
for i = 1:length(q)
    q(i) = length(unique(oa(:,i)));
end
uq = unique(q); %%设计出现的不同水平的水平数罗列
%m = histc(q,[uq-0.5;inf]); m(end) = []; %各水平因子个数罗列
m = histcounts(q,[uq-0.5;inf]); %各水平因子个数罗列
for n1 = 3:m(1)
    for n2 = 3:m(2)
        q0 = [uq(1)*ones(n1,1);uq(2)*ones(n2,1)];
        [D0,id0,aveDisc0] = MADsub_from_OA(oa,q0,'WD','out.txt');
        %导入初步计算结果
        data = importdata('out.txt');
        n0 = n1+n2;
        aveDisc = data(:,1);
        cols = data(:,2:n0+1);
        A = data(:,n0+2:2*n0+1);
        BB = data(:,2*n0+2:end);
        
        %设置精度
        k = 10; %保留到小数点后第k位
        A = round(A*10^k)/10^k;
        BB = round(BB*10^k)/10^k;
        aveDisc = round(aveDisc*10^k)/10^k;
        
        %去掉字型型相同的行
        [BB,ib] = unique(BB,'rows');
        cols = cols(ib,:);
        A = A(ib,:);
        aveDisc = aveDisc(ib,:);
        
        %按照平均均匀性排序
        [aveDisc,id] = sort(aveDisc);
        cols = cols(id,:);
        A = A(id,:);
        BB = BB(id,:);
        [~,ia] = sortrows(A);
        [~,ia] = sort(ia);
        if ia(1)~=1
            fprintf('%d %d ******\n',n1,n2);
        else
            fprintf('%d %d \n',n1,n2);
        end
        
        %结果写入文件
        %filename = strcat('Results/MA.36.3.7.6.3.','sub.',int2str(n1),'.',int2str(n2),'.txt');
        filename = strcat('Results/MA.36.2.11.3.12.','sub.',int2str(n1),'.',int2str(n2),'.txt');
        %filename = strcat('Results/MA.36.3.12.12.1.','sub.',int2str(n1),'.',int2str(n2),'.txt');
        %filename = strcat('Results/MA.72.3.24.24.1.','sub.',int2str(n1),'.',int2str(n2),'.txt');
        fid = fopen(filename,'w');
        for i = 1:n0
            fprintf(fid,'cols ');
        end
        fprintf(fid,'aveDisc Arank ');
        for i = 3:n0
            fprintf(fid,strcat('A',int2str(i)));
            fprintf(fid,' ');
        end
        for i = 1:size(BB,2)
            fprintf(fid,'BB ');
        end
        fprintf(fid,'\n');
        for i = 1:size(cols,1)
            fprintf(fid,'%d ',cols(i,:));
            fprintf(fid,'%.8f %d ',aveDisc(i),ia(i));
            fprintf(fid,'%.4f ',A(i,3:end));
            fprintf(fid,'%.4f ',BB(i,:));
            fprintf(fid,'\n');
        end
        fclose(fid);
    end
end





