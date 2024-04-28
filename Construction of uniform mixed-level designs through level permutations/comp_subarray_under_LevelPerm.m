%Ŀ�ģ����ڸ��� Mixed OA�����п����ӱ������ӱ�� aveDisc,Cols,A,BB��
%     �������в�����ͬ���ӱ�д��һ���ļ�.
%��������� Results/MA.36.2.11.3.12.sub.3.5.txt,��ʾMA.36.2.11.3.12 �����п����ӱ�
%      MA.36.2.3.3.5 �� aveDisc,Cols,A,BB �Ƚ��


%oa = importdata('../Web OA/MA.36.3.7.6.3.finney.txt');
oa = importdata('../Web OA/MA.36.2.11.3.12.txt');
%oa = importdata('../Web OA/MA.36.3.12.12.1.txt');
%oa = importdata('../Web OA/MA.72.3.24.24.1.txt');
q = zeros(size(oa,2),1);
for i = 1:length(q)
    q(i) = length(unique(oa(:,i)));
end
uq = unique(q); %%��Ƴ��ֵĲ�ͬˮƽ��ˮƽ������
%m = histc(q,[uq-0.5;inf]); m(end) = []; %��ˮƽ���Ӹ�������
m = histcounts(q,[uq-0.5;inf]); %��ˮƽ���Ӹ�������
for n1 = 3:m(1)
    for n2 = 3:m(2)
        q0 = [uq(1)*ones(n1,1);uq(2)*ones(n2,1)];
        [D0,id0,aveDisc0] = MADsub_from_OA(oa,q0,'WD','out.txt');
        %�������������
        data = importdata('out.txt');
        n0 = n1+n2;
        aveDisc = data(:,1);
        cols = data(:,2:n0+1);
        A = data(:,n0+2:2*n0+1);
        BB = data(:,2*n0+2:end);
        
        %���þ���
        k = 10; %������С������kλ
        A = round(A*10^k)/10^k;
        BB = round(BB*10^k)/10^k;
        aveDisc = round(aveDisc*10^k)/10^k;
        
        %ȥ����������ͬ����
        [BB,ib] = unique(BB,'rows');
        cols = cols(ib,:);
        A = A(ib,:);
        aveDisc = aveDisc(ib,:);
        
        %����ƽ������������
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
        
        %���д���ļ�
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





