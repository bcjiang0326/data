function enumerate_MA82441_2col(id,flag,outfilename)
%20150620 by Xiaopang
%INPUT:
%       id: MA.8.2.4.4.1.txt 选中的列
%       flag: 'CD','WD' or 'MD'
%       outfile: 输出文件的名字

%flag = 'CD'; id = [1,2,3,5];
D0 = importdata('OA from web/MA.8.2.4.4.1.txt');

D = D0(:,id);

if strcmp(flag,'CD')
    fun = @CD2_value;
elseif strcmp(flag,'WD')
    fun = @WD2_value;
elseif strcmp(flag,'MD')
    fun = @MD2_value;
else
    error('flag must be CD,WD or MD!\n');
end
n = length(id);
N = size(D0,1);
if n~=2 || id(end)~=5
    error('The length of id must be 2 and the last must be 5!\n');
end


A = perms([0,1,2,3])';
Alen = size(A,2);
B_ = Id2Design([2;2;2;2],(0:2^4-1)');
B = zeros(size(B_,1),2*size(B_,2));
for j = 1:4
    B(:,2*j-1) = B_(:,j);
    B(:,2*j) = 1-B_(:,j);
end
B = B';
Blen = size(B,2);

L = zeros(N,n);
q_lh = N*ones(n,1);

outfile = fopen(outfilename,'w');
for a1 = 1:Alen
    for a2 = 1:Alen
        k = 1;
        i1 = D(:,k)==0; i2 = D(:,k)==1;
        L(i1,k) = A(:,a1); L(i2,k) = A(:,a2)+4;
        for d1 = 1:Alen
            for d2 = 1:Blen
                k = 2;
                [~,i4] = sort(D(:,end));
                v = zeros(N,1);
                for j = 1:4
                    v(2*j-1:2*j) = 2*A(j,d1);
                end
                L(i4,k) = v+B(:,d2);
                y = fun(L,q_lh);
                fprintf(outfile,'%f\n',y);
            end
        end
    end
end
fclose(outfile);
y = importdata(outfilename);
fprintf('%d ',id); fprintf('\n');
fprintf('%.6f %.6f %.6f\n',[min(y),mean(y),std(y)]);

end
            
    








