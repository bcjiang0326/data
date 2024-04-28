function y = t_f(C,x)
n = size(x,1);
y = zeros(n,1);
for i = 1:n
    y(i) = sum(C.^x(i,:));
end
end