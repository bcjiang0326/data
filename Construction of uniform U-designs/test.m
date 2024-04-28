s = 2;
n = 8;
D0 = combins(n,s); % втсиап
N = size(D0,1);
q = s*ones(n,1);
InIter = 1000;
OutIter = 3000*InIter;
Reps = 1;
[Disc0, Disc_vec, L0] = TA_MD_LHD(D0,q,OutIter,InIter,1e-3,1e-6,Reps,0);
fprintf('%.8f\n',Disc0);
%
%fprintf('%.8f\n',MD2_value(L0,ones(n,1)*N));
%y = importdata('T_Disc_Disc0.txt');
%figure
%plot(1:size(y,1),y(:,3));
fid = fopen('256_8.txt','w');
for i = 1:N
    fprintf(fid,'%d ',L0(i,:));
    fprintf(fid,'\n');
end
fclose(fid);
%}
