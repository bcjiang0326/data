function Data_Lrg()
% 20130102 
% �� ��������� 100 ���ڵ��������ظ�������������
% ���Ƚ�����ȫ��������
datafile = fopen('LargeNoReps20121224.txt','r');
data = fscanf(datafile,'%f',[1,inf]);
fclose(datafile);

% ��������ļ�
outfile = fopen('disc.txt','w');
fprintf(outfile,'q s n & CD & q s n & CD & q s n & CD \\\\ \\hline \n');

N = length(data);
bgn = 1; %�����ݿ�ʼλ��
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