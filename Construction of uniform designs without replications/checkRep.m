function [Points, NumReps, ID] = checkRep(D,q)
% 20121210
% ��������Ƿ�����ظ������
% �㷨: 1.�������������Ӧ���ֵ����򷨵����
%       2.����Ű���������
%       3.�����Ƿ�����ظ�����Ų������ظ�����
% INPUT: 
%       D: n-by-s matrix, ����������
%       q: s-by-1 vector, design������ˮƽ��
% OUTPUT:
%       Points: t-by-s matrix, ÿ�д�������ظ��������
%       NumReps: t-by-1 vector, Points �и�������ظ�����
%       ID: t-by-1 vector, Points �и�������Ӧ���ֵ��������š�

s = size(q,2);
if s~=1
    error('CheckRep: Input variable q must be a s-by-1 vector!\n');
end
if any( q~=round(q) ) || any( q <= zeros(size(q)) )
    error('Design2Bin:The component of input variable must be positive integer!\n');
end
s = length(q);
if s~= size(D,2)
    error('CheckRep: The sizes of D and q must be matched!\n');
end
n = size(D,1);

% ��������
v = D(:,s);
temp = q(s);
for k = s-1:-1:1
    v = v + temp*(D(:,k));
    temp = temp*q(k);
end

v = sort(v);
ID = zeros(0,1);
NumReps = zeros(0,1);
i = 1;
while i < n
    if v(i) == v(i+1)
        ID = cat(1,ID,v(i));
        Num = 2;
        for j = i+2:n
            if v(j) == v(i)
                Num = Num+1;
            else
                break;
            end
        end
        NumReps = cat(1,NumReps,Num);
        i = j;
    else
        i = i+1;
    end
end

if isempty(ID)
    Points = zeros(0,s);
else
    Points = Id2Design(q,ID);
end

end

        


