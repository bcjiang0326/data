function [y,v] = minMaxValue(L,p)
% 20151015 by xiaopang
% ��������ռ� [0,1]^n �����е㵽������Ƶľ������ֵ
% ��Ƶľ��붨��Ϊ�㵽��������е����Сֵ
% L Ϊ��Ӧ��LH, size Ϊ {N,n}, Ԫ��Ϊ {0,...,N-1}
% p Ϊ���õķ��� 

epsilon = 1e-12;
[N,n] = size(L);
L = (L+0.5)/N;
V1 = combins(n,3);
V1 = (V1+0.5)/3;
y = 0;

for j = 1:size(V1,1)
    d = inf;
    %���� V1 �е� j �㵽��Ƶľ���
    for i = 1:N
        d1 = norm(V1(j,:)-L(i,:),p);
        if d1 < d-epsilon && abs(d1)>epsilon
            d = d1;
        end
    end
    if d > y+epsilon
        y = d;
        if nargout > 1
            v = V1(j,:);
        end
    end
end

V2 = sobolset(n);
V2 = net(V2,2^(16));
for j = 1:size(V2,1)
    d = inf;
    %���� V2 �е� j �㵽��Ƶľ���
    for i = 1:N
        d2 = norm(V2(j,:)-L(i,:),p);
        if d2 < d-epsilon && abs(d2)>epsilon
            d = d2;
        end
    end
    if d > y+epsilon
        y = d;
        if nargout > 1
            v = V2(j,:);
        end
    end
end
end
