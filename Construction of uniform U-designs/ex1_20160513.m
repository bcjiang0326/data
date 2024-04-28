Mlist = ['MD';'CD';'WD'];
for i = 1:3
    for n2 = 2:4
        for n1 = 2:4
            Measure = Mlist(i,:);
            fprintf('%s: \n',Measure);
            fprintf('n1=%d,  n2=%d\n',n1,n2);
            [rate_wb,rate_a,Num] = comp_ranks(n1,n2,Measure);
            fprintf('M=%d    p(r_wb)=%.3f    p(r_a)=%.3f\n',Num,rate_wb,rate_a);
            fprintf('\n');
        end
    end
end