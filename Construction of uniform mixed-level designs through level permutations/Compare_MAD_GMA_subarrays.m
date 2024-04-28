n1n2 = [3,3; 5,2; 5,3; 5,5; 4,6; 9,3; 3,8; 3,9; 3,10; 5,9];
MADcols = {[6 7 8 13 14 16];
    [5 6 7 8 9 14 16];
    [5 9 10 11 12 13 15 19];
    [2 3 6 7 9 14 15 19 20 22];
    [1 2 3 4 13 14 15 16 19 23];
    [1 2 3 4 5 6 7 8 9 13 14 15];
    [1 2 4 13 14 15 16 17 18 19 23];
    [3 5 11 13 14 16 17 18 20 21 22 23];
    [1 2 3 13 14 15 16 17 18 19 21 22 23];
    [1 4 6 7 12 13 15 16 17 18 19 20 22 23]};
GMAcols = {[5 6 9 14 16 19];
    [3 4 5 6 9 16 21];
    [1 4 5 10 11 19 20 22];
    [1 4 5 10 11 14 16 19 20 22];
    [1 4 10 11 14 16 17 19 20 22];
    [1 2 4 5 6 8 10 11 12 14 19 21];
    [1 4 9 13 14 16 17 18 20 22 23];
    [1 2 9 13 14 16 17 18 19 20 21 23];
    [1 2 7 13 15 16 17 18 19 20 21 22 23];
    [1 4 5 10 11 13 14 15 16 17 19 20 21 22]};
oa = importdata('../Web OA/MA.36.3.12.2.11.txt');

allDisc = {'WD';'CD';'MD'};


for i = 1:10
    D1 = oa(:,MADcols{i});
    A1 = GWP_asym(D1);
    Disc1 = zeros(1,3);
    T1 = zeros(1,3);
    
    D2 = oa(:,GMAcols{i});
    A2 = GWP_asym(D2);
    Disc2 = zeros(1,3);
    T2 = zeros(1,3);
    
    
    Disc1(1) = WD2_value(D1);
    Disc2(1) = WD2_value(D1);
    for j = 2:3
        flag = allDisc{j};
        T1(j) = cputime;
        Disc1(j) = minDisc_allPerms(D1,flag);
        T1(j) = cputime-T1(j);
        
        T2(j) = cputime;
        Disc2(j) = minDisc_allPerms(D2,flag);
        T2(j) = cputime-T2(j);
    end
    
    fprintf('Stage2Time(%d %d): ',n1n2(i,:));
    fprintf('%.2f  %.2f  %.2f  %.2f  %.2f  %.2f\n',T1,T2);
    %fprintf('%.6f & %.1f & %.6f & %.1f & %.6f & %.1f & %.2f \\\\\n',Disc2(1),T2(1),Disc2(2),T2(2),Disc2(3),T2(3),A2(3));
end




