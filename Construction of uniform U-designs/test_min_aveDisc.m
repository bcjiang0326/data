n1 = 3; n2 = 3;
DiscMeasure = 'MD';
D = importdata('OA from web/MA.36.3.12.2.11.txt');
D = D(:,[13:23,1:12]);


%
q0 = [2*ones(n1,1);3*ones(n2,1)];
writedata = 1;
%T0 = 1e-2; T1 = 1e-5;
T0 = 0; T1 = 0;
[aveDisc,cols,WB,A,D0] = minAveDiscSubDesign(D,q0,DiscMeasure,10,100,T0,T1,writedata);
fprintf('%.8f\n',aveDisc);
data = importdata('aveDisc.txt');
plot(1:200,data(1:200,1));
%



%[aveDisc,cols,WB,A,D0] = minAveDiscSubDesign_SA(D,q0,DiscMeasure,10,100,0,0,writedata);

%{
% ¥ ”√
minDisc = inf;
for i = 1:10^3
    cols = rand_subOA(D,q0);
    [~,aveDisc] = Weighted_wordtype_pattern(D(:,cols),q0,DiscMeasure);
    if aveDisc < minDisc-1e-11
        minDisc = aveDisc;
    end
    if mod(i,100)==0
        fprintf('%d  %.8f\n',i,minDisc);
    end
end
fprintf('%.8f\n',minDisc);
%}