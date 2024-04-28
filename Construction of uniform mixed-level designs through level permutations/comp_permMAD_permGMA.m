% 比较 Permuted MAD designs 和 Permuted GMA designs 
D = importdata('../Web OA/MA.36.2.11.3.12.txt');
for n1 = 1:11
    for n2 = 1:12
        filename = strcat('Results/MA.36.2.11.3.12.sub.',int2str(n1),'.',int2str(n2),'.txt');
        data = importdata(filename);
        j = n1+n2+2;
        if data.data(1,j)==1
            continue;
        end
        i_gma = 1;
        for i = 1:size(data.data,1)
            if data.data(i,j)==1
                i_gma = i;
                break;
            end
        end
        %计算 Permuted MAD design
        cols = data.data(1,[1:n1+n2]);
        aveCD = data.data(1,n1+n2+1);
        A3 = data.data(1,n1+n2+3);
        A4 = data.data(1,n1+n2+4);
        minCD = minCD_allPerms(D(:,cols));
        id = cols <= 11;
        cols(id) = cols(id)+12;
        cols(~id) = cols(~id)-11;
        cols = sort(cols);
        fprintf('(%d,%d)-MAD & ',[n2,n1]);
        fprintf('%d ',cols); fprintf('& ');
        fprintf('%.8f & %.8f & %.2f & %.2f \\\\ \n',[aveCD,minCD,A3,A4]);
        %计算Permuted GMA design
        cols = data.data(i_gma,[1:n1+n2]);
        aveCD = data.data(i_gma,n1+n2+1);
        A3 = data.data(i_gma,n1+n2+3);
        A4 = data.data(i_gma,n1+n2+4);
        minCD = minCD_allPerms(D(:,cols));
        id = cols <= 11;
        cols(id) = cols(id)+12;
        cols(~id) = cols(~id)-11;
        cols = sort(cols);
        fprintf('(%d,%d)-GMA & ',[n2,n1]);
        fprintf('%d ',cols); fprintf('& ');
        fprintf('%.8f & %.8f & %.2f & %.2f \\\\ \n',[aveCD,minCD,A3,A4]);
    end
end