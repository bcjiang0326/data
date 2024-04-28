%本程序检验所有 MAD subarrays 在 WD, CD 和 MD 准则下是否相同

oa = importdata('../Web OA/MA.72.4.1.6.5.3.8.2.7.txt');
allDisc = {'WD';'CD';'MD'};

fid = fopen('MA.72.GMA.MAD.distinct.txt','w');
for j = 1:3
    flag = allDisc{j};
    fprintf(fid,'%s:\n',flag);
    fprintf('%s:\n',flag);
    for n1 = 1:5
        for n2 = 1:8
            q0 = [4;6*ones(n1,1);3*ones(n2,1)];
            [~,id1,aveDisc1] = MADsub_from_OA(oa,q0,flag);
            A1 = GWP_asym(oa(:,id1),q0);
            [~,id2,A2] = GMAsub_from_OA(oa,q0);
            aveDisc2 = aveDisc_LevelPerm(oa(:,id2),q0,flag);
            if any(id1~=id2) && abs(aveDisc1-aveDisc2)>1e-10 && any(abs(A1-A2)>1e-10) 
                fprintf('(n1,n2)=(%d,%d): GMA and MAD designs are distinct!\n',n1,n2);
                fprintf(fid,'(n1,n2)=(%d,%d): id1=(',n1,n2);
                fprintf(fid,'%d ',id1);
                fprintf(fid,'), id2=(');
                fprintf(fid,'%d ',id2);
                fprintf(fid,')\n');
            else
                fprintf('(n1,n2)=(%d,%d): GMA and MAD designs are the same!\n',n1,n2);
            end
        end
    end
end
fclose(fid);

%{
WD:
(n1,n2)=(1,2): GMA and MAD designs are distinct!
(n1,n2)=(2,2): GMA and MAD designs are distinct!
(n1,n2)=(2,3): GMA and MAD designs are distinct!
(n1,n2)=(2,4): GMA and MAD designs are distinct!
(n1,n2)=(2,6): GMA and MAD designs are distinct!
(n1,n2)=(3,2): GMA and MAD designs are distinct!
(n1,n2)=(3,3): GMA and MAD designs are distinct!
(n1,n2)=(3,4): GMA and MAD designs are distinct!
(n1,n2)=(3,5): GMA and MAD designs are distinct!
(n1,n2)=(3,6): GMA and MAD designs are distinct!
(n1,n2)=(3,7): GMA and MAD designs are distinct!
(n1,n2)=(3,8): GMA and MAD designs are distinct!
(n1,n2)=(4,3): GMA and MAD designs are distinct!
(n1,n2)=(4,4): GMA and MAD designs are distinct!
(n1,n2)=(4,5): GMA and MAD designs are distinct!
(n1,n2)=(4,6): GMA and MAD designs are distinct!
(n1,n2)=(5,3): GMA and MAD designs are distinct!
(n1,n2)=(5,4): GMA and MAD designs are distinct!
(n1,n2)=(5,5): GMA and MAD designs are distinct!
(n1,n2)=(5,6): GMA and MAD designs are distinct!
(n1,n2)=(5,7): GMA and MAD designs are distinct!
CD:
(n1,n2)=(1,2): GMA and MAD designs are distinct!
(n1,n2)=(2,2): GMA and MAD designs are distinct!
(n1,n2)=(2,3): GMA and MAD designs are distinct!
(n1,n2)=(2,4): GMA and MAD designs are distinct!
(n1,n2)=(2,6): GMA and MAD designs are distinct!
(n1,n2)=(3,2): GMA and MAD designs are distinct!
(n1,n2)=(3,3): GMA and MAD designs are distinct!
(n1,n2)=(3,4): GMA and MAD designs are distinct!
(n1,n2)=(3,5): GMA and MAD designs are distinct!
(n1,n2)=(3,6): GMA and MAD designs are distinct!
(n1,n2)=(3,7): GMA and MAD designs are distinct!
(n1,n2)=(3,8): GMA and MAD designs are distinct!
(n1,n2)=(4,4): GMA and MAD designs are distinct!
(n1,n2)=(4,5): GMA and MAD designs are distinct!
(n1,n2)=(4,6): GMA and MAD designs are distinct!
(n1,n2)=(5,4): GMA and MAD designs are distinct!
(n1,n2)=(5,5): GMA and MAD designs are distinct!
(n1,n2)=(5,6): GMA and MAD designs are distinct!
(n1,n2)=(5,7): GMA and MAD designs are distinct!
MD:
(n1,n2)=(1,2): GMA and MAD designs are distinct!
(n1,n2)=(2,2): GMA and MAD designs are distinct!
(n1,n2)=(2,3): GMA and MAD designs are distinct!
(n1,n2)=(2,4): GMA and MAD designs are distinct!
(n1,n2)=(2,6): GMA and MAD designs are distinct!
(n1,n2)=(3,2): GMA and MAD designs are distinct!
(n1,n2)=(3,3): GMA and MAD designs are distinct!
(n1,n2)=(3,4): GMA and MAD designs are distinct!
(n1,n2)=(3,5): GMA and MAD designs are distinct!
(n1,n2)=(3,6): GMA and MAD designs are distinct!
(n1,n2)=(3,7): GMA and MAD designs are distinct!
(n1,n2)=(3,8): GMA and MAD designs are distinct!
(n1,n2)=(4,3): GMA and MAD designs are distinct!
(n1,n2)=(4,4): GMA and MAD designs are distinct!
(n1,n2)=(4,5): GMA and MAD designs are distinct!
(n1,n2)=(4,6): GMA and MAD designs are distinct!
(n1,n2)=(5,3): GMA and MAD designs are distinct!
(n1,n2)=(5,4): GMA and MAD designs are distinct!
(n1,n2)=(5,5): GMA and MAD designs are distinct!
(n1,n2)=(5,6): GMA and MAD designs are distinct!
(n1,n2)=(5,7): GMA and MAD designs are distinct!
%}
