function comp_betaP_PD_CD_lvl3()

N = 27 ; 
X = combins(3,3); % ������
C = [1 1 1; 1 2 0; 1 1 2; 1 0 1; 0 1 2; 1 2 2; 1 1 0]';
%D = [X,X*C];

for n = 5:10
    str1 = strcat('betaP',int2str(n),'.txt');
    str2 = strcat('PD',int2str(n),'.txt');
    str3 = strcat('CD',int2str(n),'.txt');
    outBeta = fopen(str1,'w');
    outPD = fopen(str2,'w');
    outCD = fopen(str3,'w');
    
    B = screenB(C(:,1:n-3));
    q = 3*ones(n,1);
    for k = 1:size(B,1)
        fprintf('n = %d, k = %d \n',[n,k]);
        
        D = mod( [X,X*C(:,1:n-3)+ones(N,1)*B(k,:)],3 );
        
        beta_P = beta_pattern_lvl3(D,3:10);
        CD_P = CD2_pattern(D,q,3:n);
        CD = CD2_value(D,q);
        
        fprintf(outBeta,'%.4e ',beta_P); 
        fprintf(outBeta,'%d ',[B(k,:),k] ); 
        fprintf(outBeta,'\n');
        
        fprintf(outPD,'%.6e ',CD_P); 
        fprintf(outPD,'%d ',[B(k,:),k]);
        fprintf(outPD,'\n');
    
        fprintf(outCD,'%.6e ',CD);
        fprintf(outCD,'%d ', [B(k,:),k]);
        fprintf(outCD,'\n');
    end
    fclose(outBeta); fclose(outPD); fclose(outCD);
    
    %��������
    data = importdata(str1);
    data = sortrows(data);
    outBeta = fopen(str1,'w');
    for kk=1:size(data,1)
        fprintf(outBeta,'%.4e ',data(kk,1:6) );
        fprintf(outBeta,'%d ',data(kk,7:end) );
        fprintf(outBeta,'\n');
    end
    fclose(outBeta);
    
    data = importdata(str2);
    data = sortrows(data);
    outPD = fopen(str2,'w');
    for kk=1:size(data,1)
        fprintf(outPD,'%.6e ',data(kk,1:n-2) );
        fprintf(outPD,'%d ',data(kk,n-1:end) );
        fprintf(outPD,'\n');
    end
    fclose(outPD);
    
    data = importdata(str3);
    data = sortrows(data);
    outCD = fopen(str3,'w');
    for kk=1:size(data,1)
        fprintf(outCD,'%.6e ',data(kk,1) );
        fprintf(outCD,'%d ',data(kk,2:end));
        fprintf(outCD,'\n');
    end
    fclose(outCD);
    
end

end

function B = screenB(C)
% 2013 0928
% �ο� Cheng and Ye (2004) Lemma 3.2
% Geometric isomorphism and minimum aberration for
% factorial designs with quantitative factors
% ˵����n-k �������� X, k ���������� Y = X*C+b, ���� b' ���ù��� 3^k-by-k ���� B1 
% ������Ŀ��: ���� 3^{n-k} FF ��ƣ�������(3^k+1)/2 ������ͬ����ƣ�
%       ������ɸѡ B1 �е��� b���õ����ࣨ3^k+1)/2 ���У�������� B �С�
% INPUT:
%       C: (n-k)-by-k matrix, �������� Y ���������� X ������ϵ�� 
% OUTPUT:
%       B: (3^k+1)/2-by-k matrix, ���� (3^k+1)/2 ������ͬ����ƶ�Ӧ�ĳ��� b' ���ɵľ���

k = size(C,2);
B = combins(k,3);
N = size(B,1);
for j = N:-1:2
    b = B(j,:);
    b_img = mod(2-2*sum(C)-b,3);
    for i = 1:j-1
        if all(B(i,:)==b_img)
            B = cat(1,B(1:j-1,:),B(j+1:end,:));
            break;
        end
    end
end


%{
% ���Ժ���
C= [1,1,1; 1,2,0]'; B = screenB(C);
%}
end


