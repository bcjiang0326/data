clear;clc;

oa = importdata('OA from web/MA.36.3.12.2.11.txt');
[N,n] = size(oa);
q = zeros(n,1);
for i = 1:n
    q(i) = length(unique(oa(:,i)));
end
[q,id] = sort(q);
oa = oa(:,id);

out1 = fopen('mixtest150624/ma.36.2.x.3.1.txt','a');
fprintf(out1,'measure = MD\n\n');
Meas = 'MD';
for x = 2:2
    q0 = cat(1,2*ones(x,1),3);
    %{
    [D0,A0,id0] = MA_from_asymOA(oa,q0);
    [E0,y0] = Disc_wordlength_pattern(D0,[],Meas);
    
    fprintf(out1,'MA:\n');
    fprintf(out1,'aveDisc= %.8f\n',y0);
    fprintf(out1,'id= '); fprintf(out1,'%d ',id0); fprintf(out1,'\n');
    fprintf(out1,'A= '); fprintf(out1,'%.8f ',A0); fprintf(out1,'\n');
    fprintf(out1,'E= '); fprintf(out1,'%.8f ',E0); fprintf(out1,'\n\n');
    fprintf('MA:\n');
    fprintf('aveDisc= %.8f\n',y0);
    fprintf('id= '); fprintf('%d ',id0); fprintf('\n');
    fprintf('A= '); fprintf('%.8f ',A0); fprintf('\n');
    fprintf('E= '); fprintf('%.8f ',E0); fprintf('\n\n');
    %}
    [D1,E1,id1] = DWPsub_from_OA(oa,q0,Meas);
    A1 = GWP_asym(D1);
    y1 = aveDisc_oalhd_mixlevel1(D1,q0,Meas);
    fprintf(out1,'MEA:\n');
    fprintf(out1,'aveDisc= %.8f\n',y1);
    fprintf(out1,'id= '); fprintf(out1,'%d ',id1); fprintf(out1,'\n');
    fprintf(out1,'A= '); fprintf(out1,'%.8f ',A1); fprintf(out1,'\n');
    fprintf(out1,'E= '); fprintf(out1,'%.8f ',E1); fprintf(out1,'\n\n\n');
    fprintf('MEA:\n');
    fprintf('aveDisc= %.8f\n',y1);
    fprintf('id= '); fprintf('%d ',id1); fprintf('\n');
    fprintf('A= '); fprintf('%.8f ',A1); fprintf('\n');
    fprintf('E= '); fprintf('%.8f ',E1); fprintf('\n\n\n');    
end
fclose(out1);    

%{
out2 = fopen('mixtest150624/ma.36.2.x.3.3.txt','w');
fprintf(out2,'measure = CD\n\n');
for x = 2:6
    q0 = cat(1,2*ones(x,1),[3;3;3]);
    [D0,A0,id0] = MA_from_asymOA(oa,q0);
    [E0,y0] = Disc_wordlength_pattern(D0,q0,Meas);
    
    fprintf(out2,'MA:\n');
    fprintf(out2,'aveDisc= %.8f\n',y0);
    fprintf(out2,'id= '); fprintf(out2,'%d ',id0); fprintf(out2,'\n');
    fprintf(out2,'A= '); fprintf(out2,'%.8f ',A0); fprintf(out2,'\n');
    fprintf(out2,'E= '); fprintf(out2,'%.8f ',E0); fprintf(out2,'\n\n');
    fprintf('MA:\n');
    fprintf('aveDisc= %.8f\n',y0);
    fprintf('id= '); fprintf('%d ',id0); fprintf('\n');
    fprintf('A= '); fprintf('%.8f ',A0); fprintf('\n');
    fprintf('E= '); fprintf('%.8f ',E0); fprintf('\n\n');
    
    [D1,E1,id1] = DWPsub_from_OA(oa,q0,Meas);
    A1 = GWP_asym(D1);
    y1 = aveDisc_oalhd_mixlevel1(D1,q0,Meas);
    fprintf(out2,'MEA:\n');
    fprintf(out2,'aveDisc= %.8f\n',y1);
    fprintf(out2,'id= '); fprintf(out2,'%d ',id1); fprintf(out2,'\n');
    fprintf(out2,'A= '); fprintf(out2,'%.8f ',A1); fprintf(out2,'\n');
    fprintf(out2,'E= '); fprintf(out2,'%.8f ',E1); fprintf(out2,'\n\n\n'); 
    fprintf('MEA:\n');
    fprintf('aveDisc= %.8f\n',y1);
    fprintf('id= '); fprintf('%d ',id1); fprintf('\n');
    fprintf('A= '); fprintf('%.8f ',A1); fprintf('\n');
    fprintf('E= '); fprintf('%.8f ',E1); fprintf('\n\n\n');
end
fclose(out2);
%}




