function [DiscMin,DiscVec,Dopt] = minDisc_allPerms(D0,flag)
%20181128 by Bochuan Jiang
%�������һ��balance design��ȫ��ˮƽ�û�����С��CDֵ
%Input:
%   D: balance design
%   flag: 'CD','WD' or 'MD'
%Output:
%   DiscMin: ��С discrepancy.
%   DiscVec: ���� permuted D �� discrepancy.
%   Dopt: ���ŵ� permuted D.
if ~isBalanced(D0)
    error('D is not balance!\n');
end
if ~strcmp(flag,'CD') && ~strcmp(flag,'WD') && ~strcmp(flag,'MD')
    error('Wrong flag!\n');
end
[~,n] = size(D0);
q0 = zeros(n,1);
for j = 1:n
    q0(j) = length(unique(D0(:,j)));
end
allPerms = cell(n,1);
allPerms_size = zeros(n,1);
for j = 1:n
    if q0(j) == 2
        allPerms{j} = [0,1];
    elseif q0(j) == 3
        allPerms{j} = [0,1,2;1,2,0;2,0,1];
    else
        allPerms{j} = sortrows(perms(0:q0(j)-1));
        Nj = size(allPerms{j},1);
        v = zeros(1,Nj/2);
        i3 = 1;
        for i1 = 1:Nj-1
            for i2 = i1+1:Nj
                if all(allPerms{j}(i2,:)==fliplr(allPerms{j}(i1,:)))
                    v(i3) = i2;
                    i3 = i3+1;
                    break;
                end
            end
        end
        allPerms{j}(v,:) = [];                
    end
    allPerms_size(j) = size(allPerms{j},1);
end
D = D0;
DiscMin = inf;
DiscVec = zeros(0,1);
if nargout >= 3
    Dopt = D;
end
isdisp = 1; %�Ƿ�����׶ν�������û�ָ��
if isdisp
    t0 = cputime;
end
for k = 0:prod(allPerms_size)-1
    vec = Id2Design(allPerms_size,k)+1;
    for j = 1:n %�û���j��
        for L = 0:q0(j)-1 %�û�ˮƽL
            D(D0(:,j)==L,j) = allPerms{j}(vec(j),L+1);
        end
    end
    DiscVec = cat(1,DiscVec,Disc2_value(D,q0,flag));
    if DiscVec(end) < DiscMin-1e-12
        DiscMin = DiscVec(end);
        if nargout >= 3
            Dopt = D;
        end
    end
    
    %���ڳ���1e5�����Ҫ��ʾ����
    disp_j = n-1;
    prod_size = allPerms_size(n);
    while prod_size < 1e5 && disp_j > 0
        prod_size = prod_size*allPerms_size(disp_j);
        disp_j = disp_j-1;
    end
    if prod_size > 1e5 && disp_j == 0
        disp_j = 1;
    end
    if disp_j > 0
        if all(vec(disp_j+1:end)==1)
            t1 = cputime;
            fprintf('%d ',vec(1:disp_j));
            fprintf(' time = %.2f\n',t1-t0);
            t0 = t1;
        end
    end
end
end
%{
D = importdata('../Web OA/MA.36.3.7.6.3.finney.txt');
D1 = D(:,[1,2,4,5,10]);
D2 = D(:,[1,3,5,7,8]);
[CDmin1,CDvec1,Dopt1] = minDisc_allPerms(D1,'CD');
[CDmin2,CDvec2,Dopt2] = minDisc_allPerms(D2,'CD');
%}
