function [c1,c2,a1,a2,cc1,cc2,aa1,aa2] = WDcoef(N,s)
b = zeros(s,s);
bb = zeros(s,s);
a = N/s;
for l1 = 1:s
    for l2 = 1:s
        for k1 = (l1-1)*a:l1*a-1
            for k2 = (l2-1)*a:l2*a-1
                if k1~=k2
                    b(l1,l2) = b(l1,l2)+f(N,k1,k2);
                    bb(l1,l2) = bb(l1,l2)+ff(N,k1,k2);
                end
            end
        end
    end
end
c1 = 0;c2 = 0;
cc1 = 0;cc2 = 0;
for i = 1:s
    c2 = c2+b(i,i);
    cc2 = cc2+bb(i,i);
    for j = 1:s
        if i~=j
            c1 = c1+b(i,j);
            cc1 = cc1+bb(i,j);
        end
    end
end
c2 = c2*a*(s-1)/(a-1);
cc2 = cc2*a*(s-1)/(a-1);
a1 = (c2+c1*(s-1))/N^2/(s-1);
a2 = (c2-c1)/(c2+c1*(s-1));
aa1 = (cc2+cc1*(s-1))/N^2/(s-1);
aa2 = (cc2-cc1)/(cc2+cc1*(s-1));
end

function y = f(N,i,j)
y = 1.5-abs(i-j)/N+(i-j)^2/N^2;
end

function y = ff(N,i,j)
y = f(N,i,j)-1;
end