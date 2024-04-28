function search_combines_betaPtn(n,m,s)
% 2013 0929
% ������n ���������� i ������ˮƽΪ s_i�������ṻ0,1, ... , s_i-1 �����ɶ�
% ������Ҫ m �����ɶȣ�������п��ܵķ�����
% ���̽��ܣ���������õݹ�ķ���������ȼ��������
% INPUT:
%           
%           n: ��ǰ��ѡ��ı�������
%           m: ��ǰ�������ɶȵĸ���
%           s: scalar or n-by-1 vector, ��Ǹ���������ˮƽ
if n~=round(n) || n < 1
    error('n must be a positive integer!\n');
end
if length(s)==1
    s = s*ones(1,n);
elseif length(s)~=n
    error('s and n do not match!\n');
end
if any(s~=round(s)) || m~=round(m) || min(s)<1 || m<1
    error(' The components of all inputs must be positive integers!\n');
end

b = zeros(1,n);
rec_Ntree(b,n,m,s);

end

function rec_Ntree(b,n,m,s)
%INPUT:   
%   b: 1-by-n vector, ��Ǹ��������ṩ�����ɶ�
if sum(s)-n<m
    return;
end
if m == 0
    fprintf('%d ',b);
    fprintf('\n');
    return;
end

for k = s(n)-1:-1:0
    if m >= k
        b(n) = k;
        rec_Ntree(b,n-1,m-k,s(1:n-1));
    end
end
end

%{
% ���Ժ���1
n = 5; m = 3; s = 3;
search_combines(n,m,s);
% ���Ժ���2
n = 5; m = 3; s = [3,2,3,2,3];
search_combines(n,m,s);
%}


