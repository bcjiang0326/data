flag = 'WD';
s = 3;
oa = importdata('OA from web/oa.27.13.3.2.txt');

%MA.27.3.12
D0 = oa(:,1:end-1); 
[N,n] = size(D0);
WD1 = WD2_value(D0);
projWD1 = Proj_Disc2(2,D0,flag);
A1 = GWP(D0);
UD = importdata(strcat('UD from web/WD/Level3/',int2str(n),'_',int2str(N),'.txt'));
WD2 = WD2_value(UD);
projWD2 = Proj_Disc2(2,UD,flag);
A2 = GWP(UD);
y = [N,s,n,WD1,projWD1,A1(2),WD2,projWD2,A2(2)];
fprintf('%d & %d & %d & & %.2f & %.6f & %.2f & & %.2f & %.6f & %.2f \n',y);
%a = (16*s^2+1)/36/s^4+2/n/(n-1)*((s+1)/6/s^2)^2*A2(2);
%fprintf('%.6f\n',a);

%MA.27.3.13
D1 = oa;
[N,n] = size(D1);
WD1 = WD2_value(D1);
projWD1 = Proj_Disc2(2,D1,flag);
A1 = GWP(D1);
UD = importdata(strcat('UD from web/WD/Level3/',int2str(n),'_',int2str(N),'.txt'));
UD = UD-1;
WD2 = WD2_value(UD);
projWD2 = Proj_Disc2(2,UD,flag);
A2 = GWP(UD);
y = [N,s,n,WD1,projWD1,A1(2),WD2,projWD2,A2(2)];
fprintf('%d & %d & %d & & %.2f & %.6f & %.2f & & %.2f & %.6f & %.2f \n',y);

%MA.27.3.14
D2 = [oa,oa(:,1)];
[N,n] = size(D2);
WD1 = WD2_value(D2);
projWD1 = Proj_Disc2(2,D2,flag);
A1 = GWP(D2);
UD = importdata(strcat('UD from web/WD/Level3/',int2str(n),'_',int2str(N),'.txt'));
UD = UD-1;
WD2 = WD2_value(UD);
projWD2 = Proj_Disc2(2,UD,flag);
A2 = GWP(UD);
y = [N,s,n,WD1,projWD1,A1(2),WD2,projWD2,A2(2)];
fprintf('%d & %d & %d & & %.2f & %.6f & %.2f & & %.2f & %.6f & %.2f \n',y);

%MA.18.3.12
D3 = oa(oa(:,end)~=2,1:end-1);
[N,n] = size(D3);
WD1 = WD2_value(D3);
projWD1 = Proj_Disc2(2,D3,flag);
A1 = GWP(D3);
UD = importdata(strcat('UD from web/WD/Level3/',int2str(n),'_',int2str(N),'.txt'));
UD = UD-1;
WD2 = WD2_value(UD);
projWD2 = Proj_Disc2(2,UD,flag);
A2 = GWP(UD);
y = [N,s,n,WD1,projWD1,A1(2),WD2,projWD2,A2(2)];
fprintf('%d & %d & %d & & %.2f & %.6f & %.2f & & %.2f & %.6f & %.2f \n',y);

%MA.18.3.8
D4 = importdata('MAD/level3/MAD.18.3.8.txt');
[N,n] = size(D4);
WD1 = WD2_value(D4);
projWD1 = Proj_Disc2(2,D4,flag);
A1 = GWP(D4);
UD = importdata(strcat('UD from web/WD/Level3/',int2str(n),'_',int2str(N),'.txt'));
UD = UD-1;
WD2 = WD2_value(UD);
projWD2 = Proj_Disc2(2,UD,flag);
A2 = GWP(UD);
y = [N,s,n,WD1,projWD1,A1(2),WD2,projWD2,A2(2)];
fprintf('%d & %d & %d & & %.2f & %.6f & %.2f & & %.2f & %.6f & %.2f \n',y);

%MA.18.3.16
D5 = importdata('MAD/level3/MAD.18.3.16.txt');
[N,n] = size(D5);
WD1 = WD2_value(D5);
projWD1 = Proj_Disc2(2,D5,flag);
A1 = GWP(D5);
UD = importdata(strcat('UD from web/WD/Level3/',int2str(n),'_',int2str(N),'.txt'));
UD = UD-1;
WD2 = WD2_value(UD);
projWD2 = Proj_Disc2(2,UD,flag);
A2 = GWP(UD);
y = [N,s,n,WD1,projWD1,A1(2),WD2,projWD2,A2(2)];
fprintf('%d & %d & %d & & %.2f & %.6f & %.2f & & %.2f & %.6f & %.2f \n',y);

%MA.18.3.17
D6 = importdata('MAD/level3/MAD.18.3.17.txt');
[N,n] = size(D6);
WD1 = WD2_value(D6);
projWD1 = Proj_Disc2(2,D6,flag);
A1 = GWP(D6);
UD = importdata(strcat('UD from web/WD/Level3/',int2str(n),'_',int2str(N),'.txt'));
UD = UD-1;
WD2 = WD2_value(UD);
projWD2 = Proj_Disc2(2,UD,flag);
A2 = GWP(UD);
y = [N,s,n,WD1,projWD1,A1(2),WD2,projWD2,A2(2)];
fprintf('%d & %d & %d & & %.2f & %.6f & %.2f & & %.2f & %.6f & %.2f \n',y);

%MA.18.3.18
D7 = importdata('MAD/level3/MAD.18.3.18.txt');
[N,n] = size(D7);
WD1 = WD2_value(D7);
projWD1 = Proj_Disc2(2,D7,flag);
A1 = GWP(D7);
UD = importdata(strcat('UD from web/WD/Level3/',int2str(n),'_',int2str(N),'.txt'));
UD = UD-1;
WD2 = WD2_value(UD);
projWD2 = Proj_Disc2(2,UD,flag);
A2 = GWP(UD);
y = [N,s,n,WD1,projWD1,A1(2),WD2,projWD2,A2(2)];
fprintf('%d & %d & %d & & %.2f & %.6f & %.2f & & %.2f & %.6f & %.2f \n',y);