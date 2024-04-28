%20151015 Xiaopang
% 比较 OA(36,2^11 3^12,2) 下，TA 和 MA-based 算法差异 
clear;clc;
epsilon = 1e-10; T0 = 1e-2;  T1 = 1e-6; 
InIter = 1e4; OutIter = 100*InIter; Reps = 30;
%InIter = 1e3; OutIter = 10*InIter; Reps = 3;


OA = importdata('OA from web/MA.36.3.12.2.11.txt');
[N,ncols] = size(OA);
q = zeros(ncols,1);
for i = 1:ncols
    q(i) = length(unique(OA(:,i)));
end
[q,id] = sort(q);
OA = OA(:,id);
%n1 = 11; n2 = 12;


outfile = fopen('MA LHD 20151024 MD oa36 mwa ta/1oa36_20151024.txt','w');
fprintf(outfile,'ma.36.2.11.3.12.txt');
fprintf(outfile,'   &   &   MWA-based &   TA  \\\\ \n');
fprintf(outfile,'n_1 & n_2 & AMD & time & AMD & time\\\\ \n');

%21%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q_oa =[2;2;3];
id=[1 2 12];
D = OA(:,id);
n10 = sum(q_oa==2);
n20 = sum(q_oa==3);
n = length(q_oa);
q_lh = N*ones(n,1);
Disc_vec0 = zeros(Reps,1);
t0 = cputime;
for rep = 1:Reps
    [Disc_vec0(rep), ~, L0] = TA_MD_LHD(D,q_oa,OutIter,InIter,T0,T1,1,0,1);
    Lname = strcat('MA LHD 20151024 MD oa36 mwa ta/L_mwa',int2str(n10),int2str(n20),'_',int2str(rep));
    out= fopen(Lname,'w');
    for i = 1:N
        fprintf(out,'%d ',L0(i,:));
        fprintf(out,'\n');
    end
    fclose(out);
end
t0 = cputime-t0;

y = [n10,n20,mean(Disc_vec0),t0, nan, nan];
fprintf(outfile,'%d   %d   %.4e    %.1f   %.4e    %.1f  \n',y);
fprintf('%d  &  %d  &  %.4e &  %.1f   &  %.4e &  %.1f  \\\\ \n',y);


%12%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q_oa =[2;3;3];
id=[1 12 13];
D = OA(:,id);
n10 = sum(q_oa==2);
n20 = sum(q_oa==3);
n = length(q_oa);
q_lh = N*ones(n,1);
t0 = cputime;
for rep = 1:Reps
    [Disc_vec0(rep), ~, L0] = TA_MD_LHD(D,q_oa,OutIter,InIter,T0,T1,1,0,1);
    Lname = strcat('MA LHD 20151024 MD oa36 mwa ta/L_mwa',int2str(n10),int2str(n20),'_',int2str(rep));
    out= fopen(Lname,'w');
    for i = 1:N
        fprintf(out,'%d ',L0(i,:));
        fprintf(out,'\n');
    end
    fclose(out);
end
t0 = cputime-t0;

% 采用 TA 算法 计算 uniform LHD
t2 = cputime;
Disc_vec2 = zeros(Reps,1);
for rep = 1:Reps
    [Disc_vec2(rep),~, L2] = TA_MD(N,q_lh,OutIter*1.4,InIter*1.4,T0,T1,1,0);
    Lname = strcat('MA LHD 20151024 MD oa36 mwa ta/L_ta_',int2str(n),'_',int2str(rep));
    out= fopen(Lname,'w');
    for i = 1:N
        fprintf(out,'%d ',L2(i,:));
        fprintf(out,'\n');
    end
    fclose(out);
end
t2 = cputime-t2;

y = [n10,n20,mean(Disc_vec0),t0,mean(Disc_vec2),t2];
fprintf(outfile,'%d  %d    %.4e    %.1f    %.4e   %.1f   \n',y);
fprintf('%d  &  %d  &  %.4e &  %.1f  &  %.4e &  %.1f   \\\\ \n',y);

%22%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q_oa =[2;2;3;3];
id=[1 4 12 13];
D = OA(:,id);
n10 = sum(q_oa==2);
n20 = sum(q_oa==3);
n = length(q_oa);
q_lh = N*ones(n,1);
t0 = cputime;
Disc_vec0 = zeros(Reps,1);
for rep = 1:Reps
    [Disc_vec0(rep), ~, L0] = TA_MD_LHD(D,q_oa,OutIter,InIter,T0,T1,1,0,1);
    Lname = strcat('MA LHD 20151024 MD oa36 mwa ta/L_mwa',int2str(n10),int2str(n20),'_',int2str(rep));
    out= fopen(Lname,'w');
    for i = 1:N
        fprintf(out,'%d ',L0(i,:));
        fprintf(out,'\n');
    end
    fclose(out);
end
t0 = cputime-t0;


    
% 3.直接采用 TA 算法 计算 uniform LHD
t2 = cputime;
Disc_vec2 = zeros(Reps,1);
for rep = 1:Reps
    [Disc_vec2(rep),~, L2] = TA_MD(N,q_lh,OutIter*1.4,InIter*1.4,T0,T1,1,0);
    Lname = strcat('MA LHD 20151024 MD oa36 mwa ta/L_ta_',int2str(n),'_',int2str(rep));
    out= fopen(Lname,'w');
    for i = 1:N
        fprintf(out,'%d ',L2(i,:));
        fprintf(out,'\n');
    end
    fclose(out);
end
t2 = cputime-t2;

y = [n10,n20,mean(Disc_vec0),t0,mean(Disc_vec2),t2];
fprintf(outfile,'%d  %d    %.4e    %.1f    %.4e   %.1f   \n',y);
fprintf('%d  &  %d  &  %.4e &  %.1f  &  %.4e &  %.1f   \\\\ \n',y);


%32%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q_oa =[2;2;2;3;3];
id=[1 2 7 12 16 ];
D = OA(:,id);
n10 = sum(q_oa==2);
n20 = sum(q_oa==3);
n = length(q_oa);
q_lh = N*ones(n,1);
t0 = cputime;
for rep = 1:Reps
    [Disc_vec0(rep), ~, L0] = TA_MD_LHD(D,q_oa,OutIter,InIter,T0,T1,1,0,1);
    Lname = strcat('MA LHD 20151024 MD oa36 mwa ta/L_mwa',int2str(n10),int2str(n20),'_',int2str(rep));
    out= fopen(Lname,'w');
    for i = 1:N
        fprintf(out,'%d ',L0(i,:));
        fprintf(out,'\n');
    end
    fclose(out);
end
t0 = cputime-t0;

y = [n10,n20,mean(Disc_vec0),t0, nan, nan];
fprintf(outfile,'%d   %d   %.4e    %.1f   %.4e    %.1f  \n',y);
fprintf('%d  &  %d  &  %.4e &  %.1f   &  %.4e &  %.1f  \\\\ \n',y);

%23%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q_oa =[2;2;3;3;3];
id=[2 8 12 15 20 ];
D = OA(:,id);
n10 = sum(q_oa==2);
n20 = sum(q_oa==3);
n = length(q_oa);
q_lh = N*ones(n,1);
t0 = cputime;
for rep = 1:Reps
    [Disc_vec0(rep), ~, L0] = TA_MD_LHD(D,q_oa,OutIter,InIter,T0,T1,1,0,1);
    Lname = strcat('MA LHD 20151024 MD oa36 mwa ta/L_mwa',int2str(n10),int2str(n20),'_',int2str(rep));
    out= fopen(Lname,'w');
    for i = 1:N
        fprintf(out,'%d ',L0(i,:));
        fprintf(out,'\n');
    end
    fclose(out);
end
t0 = cputime-t0;

% 采用 TA 算法 计算 uniform LHD
t2 = cputime;
Disc_vec2 = zeros(Reps,1);
for rep = 1:Reps
    [Disc_vec2(rep),~, L2] = TA_MD(N,q_lh,OutIter*1.4,InIter*1.4,T0,T1,1,0);
    Lname = strcat('MA LHD 20151024 MD oa36 mwa ta/L_ta_',int2str(n),'_',int2str(rep));
    out= fopen(Lname,'w');
    for i = 1:N
        fprintf(out,'%d ',L2(i,:));
        fprintf(out,'\n');
    end
    fclose(out);
end
t2 = cputime-t2;

y = [n10,n20,mean(Disc_vec0),t0,mean(Disc_vec2),t2];
fprintf(outfile,'%d  %d    %.4e    %.1f    %.4e   %.1f   \n',y);
fprintf('%d  &  %d  &  %.4e &  %.1f  &  %.4e &  %.1f   \\\\ \n',y);

%33%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q_oa = [2;2;2;3;3;3];
id = [1 2 4 17 18 19];
D = OA(:,id);
n10 = sum(q_oa==2);
n20 = sum(q_oa==3);
n = length(q_oa);
t0 = cputime;
for rep = 1:Reps
    [Disc_vec0(rep), ~, L0] = TA_MD_LHD(D,q_oa,OutIter,InIter,T0,T1,1,0,1);
    Lname = strcat('MA LHD 20151024 MD oa36 mwa ta/L_mwa',int2str(n10),int2str(n20),'_',int2str(rep));
    out= fopen(Lname,'w');
    for i = 1:N
        fprintf(out,'%d ',L0(i,:));
        fprintf(out,'\n');
    end
    fclose(out);
end
t0 = cputime-t0;

% 采用 TA 算法 计算 uniform LHD
n10 = 3; 
n20 = 3; 
n = n10+n20;
q_lh = N*ones(n,1);
t2 = cputime;
Disc_vec2 = zeros(Reps,1);
for rep = 1:Reps
    [Disc_vec2(rep),~, L2] = TA_MD(N,q_lh,OutIter*1.4,InIter*1.4,T0,T1,1,0);
    Lname = strcat('MA LHD 20151024 MD oa36 mwa ta/L_ta_',int2str(n),'_',int2str(rep));
    out= fopen(Lname,'w');
    for i = 1:N
        fprintf(out,'%d ',L2(i,:));
        fprintf(out,'\n');
    end
    fclose(out);
end
t2 = cputime-t2;

y = [n10,n20,mean(Disc_vec0),t0, mean(Disc_vec2),t2];
fprintf(outfile,'%d  %d    %.4e    %.1f    %.4e   %.1f   \n',y);
fprintf('%d  &  %d  &  %.4e &  %.1f  &  %.4e &  %.1f   \\\\ \n',y);

%43%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q_oa =[2;2;2;2;3;3;3];
id=[2 4 7 9 16 17 20 ];
D = OA(:,id);
n10 = sum(q_oa==2);
n20 = sum(q_oa==3);
n = length(q_oa);
q_lh = N*ones(n,1);
t0 = cputime;
for rep = 1:Reps
    [Disc_vec0(rep), ~, L0] = TA_MD_LHD(D,q_oa,OutIter,InIter,T0,T1,1,0,1);
    Lname = strcat('MA LHD 20151024 MD oa36 mwa ta/L_mwa',int2str(n10),int2str(n20),'_',int2str(rep));
    out= fopen(Lname,'w');
    for i = 1:N
        fprintf(out,'%d ',L0(i,:));
        fprintf(out,'\n');
    end
    fclose(out);
end
t0 = cputime-t0;

y = [n10,n20,mean(Disc_vec0),t0, nan, nan];
fprintf(outfile,'%d   %d   %.4e    %.1f   %.4e    %.1f  \n',y);
fprintf('%d  &  %d  &  %.4e &  %.1f   &  %.4e &  %.1f  \\\\ \n',y);


%24%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q_oa =[2;2;3;3;3;3];
id = [2 8 12 15 21 22]; 
D = OA(:,id);
n10 = sum(q_oa==2);
n20 = sum(q_oa==3);
n = length(q_oa);
q_lh = N*ones(n,1);
Disc_vec0 = zeros(Reps,1);
t0 = cputime;
for rep = 1:Reps
    [Disc_vec0(rep), ~, L0] = TA_MD_LHD(D,q_oa,OutIter,InIter,T0,T1,1,0,1);
    Lname = strcat('MA LHD 20151024 MD oa36 mwa ta/L_mwa',int2str(n10),int2str(n20),'_',int2str(rep));
    out= fopen(Lname,'w');
    for i = 1:N
        fprintf(out,'%d ',L0(i,:));
        fprintf(out,'\n');
    end
    fclose(out);
end
t0 = cputime-t0;

y = [n10,n20,mean(Disc_vec0),t0, nan, nan];
fprintf(outfile,'%d   %d   %.4e    %.1f   %.4e    %.1f  \n',y);
fprintf('%d  &  %d  &  %.4e &  %.1f   &  %.4e &  %.1f  \\\\ \n',y);


%34%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q_oa =[2;2;2;3;3;3;3];
id=[1 2 3 16 21 22 23];
D = OA(:,id);
n10 = sum(q_oa==2);
n20 = sum(q_oa==3);
n = length(q_oa);
q_lh = N*ones(n,1);
t0 = cputime;
for rep = 1:Reps
    [Disc_vec0(rep), ~, L0] = TA_MD_LHD(D,q_oa,OutIter,InIter,T0,T1,1,0,1);
    Lname = strcat('MA LHD 20151024 MD oa36 mwa ta/L_mwa',int2str(n10),int2str(n20),'_',int2str(rep));
    out= fopen(Lname,'w');
    for i = 1:N
        fprintf(out,'%d ',L0(i,:));
        fprintf(out,'\n');
    end
    fclose(out);
end
t0 = cputime-t0;

% 采用 TA 算法 计算 uniform LHD
t2 = cputime;
Disc_vec2 = zeros(Reps,1);
for rep = 1:Reps
    [Disc_vec2(rep),~, L2] = TA_MD(N,q_lh,OutIter*1.4,InIter*1.4,T0,T1,1,0);
    Lname = strcat('MA LHD 20151024 MD oa36 mwa ta/L_ta_',int2str(n),'_',int2str(rep));
    out= fopen(Lname,'w');
    for i = 1:N
        fprintf(out,'%d ',L2(i,:));
        fprintf(out,'\n');
    end
    fclose(out);
end
t2 = cputime-t2;

y = [n10,n20,mean(Disc_vec0),t0,mean(Disc_vec2),t2];
fprintf(outfile,'%d  %d    %.4e    %.1f    %.4e   %.1f   \n',y);
fprintf('%d  &  %d  &  %.4e &  %.1f  &  %.4e &  %.1f   \\\\ \n',y);


%44%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
id=[1 2 3 7 16 21 22 23];
q_oa =[2;2;2;2;3;3;3;3];
D = OA(:,id);
n10 = sum(q_oa==2);
n20 = sum(q_oa==3);
n = length(q_oa);
q_lh = N*ones(n,1);
t0 = cputime;
for rep = 1:Reps
    [Disc_vec0(rep), ~, L0] = TA_MD_LHD(D,q_oa,OutIter,InIter,T0,T1,1,0,1);
    Lname = strcat('MA LHD 20151024 MD oa36 mwa ta/L_mwa',int2str(n10),int2str(n20),'_',int2str(rep));
    out= fopen(Lname,'w');
    for i = 1:N
        fprintf(out,'%d ',L0(i,:));
        fprintf(out,'\n');
    end
    fclose(out);
end
t0 = cputime-t0;

% 采用 TA 算法 计算 uniform LHD
t2 = cputime;
Disc_vec2 = zeros(Reps,1);
for rep = 1:Reps
    [Disc_vec2(rep),~, L2] = TA_MD(N,q_lh,OutIter*1.4,InIter*1.4,T0,T1,1,0);
    Lname = strcat('MA LHD 20151024 MD oa36 mwa ta/L_ta_',int2str(n),'_',int2str(rep));
    out= fopen(Lname,'w');
    for i = 1:N
        fprintf(out,'%d ',L2(i,:));
        fprintf(out,'\n');
    end
    fclose(out);
end
t2 = cputime-t2;

y = [n10,n20,mean(Disc_vec0),t0,mean(Disc_vec2),t2];
fprintf(outfile,'%d  %d    %.4e    %.1f    %.4e   %.1f   \n',y);
fprintf('%d  &  %d  &  %.4e &  %.1f  &  %.4e &  %.1f   \\\\ \n',y);




