%本程序检验所有 MAD subarrays 在 WD, CD 和 MD 准则下是否相同

oa = importdata('../Web OA/MA.72.4.1.6.5.3.8.2.7.txt');
allDisc = {'WD';'CD';'MD'};

for n1 = 1:5
    for n2 = 1:8
        q0 = [4;6*ones(n1,1);3*ones(n2,1)];
        T0 = zeros(1,3);
        T0(1) = cputime;
        [~,id1] = MADsub_from_OA(oa,q0,'WD');
        T0(1) = cputime-T0(1);
        T0(2) = cputime;
        [~,id11] = MADsub_from_OA(oa,q0,'CD');
        T0(2) = cputime-T0(2);
        T0(3) = cputime;
        [~,id111] = MADsub_from_OA(oa,q0,'MD');
        T0(3) = cputime-T0(3);
        
        if any(id1~=id11) || any(id1~=id111) || any(id11~=id111)
            fprintf('(n1,n2)=(%d,%d): MAD designs are distinct!\n',n1,n2);
        else
            fprintf('(n1,n2)=(%d,%d): MAD designs are the same!\n',n1,n2);
        end        
    end
end

%{
(n1,n2)=(4,3): MAD designs are distinct!
(n1,n2)=(5,3): MAD designs are distinct!
(n1,n2)=(5,4): MAD designs are distinct!
%}
