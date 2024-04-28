function AverageLoops(levels,s)

k = 1;
while 0.5*levels^k<60
    k = k+1;
end

zone = levels:levels^(s-k):0.5*levels^s;


Iter = 100000;

q = levels*ones(s,1);
N = levels^s;
dd = ones(s,1);
for i=1:s-1
    dd(i)=prod(q(i+1:end));
end

fname = strcat('ave_loops_q',int2str(levels),'s',int2str(s),'.txt');
outfile0 = fopen(fname,'w');
for n = zone
    outfile = fopen('rep.txt','w');
    
    [D,ID] = rand_U_Type_orth_design(q,n);
    for Oiter = 1:Iter
        %����Ҫ�������м���������Ԫ��    
        isrepeat = 0;
        k=randi(s);
        i=randi(n);
        j=randi(n);
        while j==i || D(i,k)==D(j,k)    
            j=randi(n);
        end
        if i > j % ǿ�� i Ϊ��С���Ǹ������淽��
            temp=i; i=j; j=temp;
        end
        IDi = ID(i) + (D(j,k)-D(i,k))*dd(k); %�˴���¼�´˶�����
        IDj = ID(j) + (D(i,k)-D(j,k))*dd(k); %���������� ID
        if ~isempty(find(ID==IDi,1)) || ~isempty(find(ID==IDj,1))        
            isrepeat = 1;        
        end        
        %��ʱ���ȿ��ǽ����Ƿ������ظ���
        %�������ظ��У������½��н���
        rep = 0;
        while isrepeat
            rep = rep +1;
            isrepeat = 0;
            k=randi(s);
            i=randi(n);
            j=randi(n);
            while j==i || D(i,k)==D(j,k)
                j=randi(n);
            end
            if i > j
                temp=i; i=j; j=temp;
            end
            %�����Ƿ���ܳ����ظ�
            IDi = ID(i) + (D(j,k)-D(i,k))*dd(k);
            IDj = ID(j) + (D(i,k)-D(j,k))*dd(k);
            if ~isempty(find(ID==IDi,1)) || ~isempty(find(ID==IDj,1))
                isrepeat = 1;
            end
        end
        % ���� D,ID
        temp=D(i,k);D(i,k)=D(j,k);D(j,k)=temp;
        ID(i) = IDi; ID(j) = IDj;
        
        fprintf(outfile,'%d ',rep);
    end
    fclose(outfile);
    rep = importdata('rep.txt');
    u = [n,mean(rep),max(rep),std(rep)];
    cout_8f([n,mean(rep),max(rep),std(rep)]);
    fprintf(outfile0,'%d %.8f %d %.8f\n',u);
end
fclose(outfile0);

end

%{
% ���Դ���
levels = 2; s = 8;  AverageLoops(levels,s);
q = levels*ones(s,1); N = levels^s;
fname = strcat('ave_loops_q',int2str(levels),'s',int2str(s),'.txt');
result = importdata(fname);
plot(result(:,1)/N,result(:,2),'-','LineWidth',0.01); hold on;

levels = 3; s = 7;  AverageLoops(levels,s);
q = levels*ones(s,1); N = levels^s;
fname = strcat('ave_loops_q',int2str(levels),'s',int2str(s),'.txt');
result = importdata(fname);
plot(result(:,1)/N,result(:,2),'-.','LineWidth',0.01);hold on;

levels = 5; s = 5;  AverageLoops(levels,s);
q = levels*ones(s,1); N = levels^s;
fname = strcat('ave_loops_q',int2str(levels),'s',int2str(s),'.txt');
result = importdata(fname);
plot(result(:,1)/N,result(:,2),'--','LineWidth',0.01); hold on;

set(gca,'YGrid','on');
legend('q=2,m=8','q=3,m=7','q=5,m=5',2);
xlabel('n/N');
ylabel('Average number of loops');
%}
    