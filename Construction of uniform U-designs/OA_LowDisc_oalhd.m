function [oa,A0,PD,Disc2,id0] = OA_LowDisc_oalhd(OA,n,disc)
%20150603 by Xiaopang
%在现有的OA中搜索具有n列的sub OA，使得其生成的OA-based LH 具有最小的average squared discrepancy
%Input:
%   OA: 正交表，猎术大于n
%   n: 待选区的列数
%   disc: string 型, 'CD','WD','MD'
%Output:
%   oa: sub OA
%   A0: GWP of sub OA, i.e., (A_1(oa),...,A_n(oa))
%   PD: (PD_1,...,PD_n), projection discrepancy pattern 参见Hickernell and Liu (2002)
%   Disc2: CD2,WD2 or MD2 values of sub OA
%   id0: OA(:,id0) = oa;

[N,n1] = size(OA);
if n>n1
    error('Input error!\n');
end
epsilon = 1e-12;
s = length(unique(OA(:,1)));
if strcmp(disc,'CD')
    [h,b,H,B] = CD_GWPoa_lhd(N,n,s);
elseif strcmp(disc,'WD')
    [h,b,H,B] = WD_GWPoa_lhd(N,n,s);
elseif strcmp(disc,'MD')
    [h,b,H,B] = MD_GWPoa_lhd(N,n,s);
else
    error('disc must be CD, WD or MD!\n');
end


id = 1:n;
P = Krawtchouk_polyMat(n,s);
A0 = inf*ones(1,n);
Disc2 = inf;
while id(1)~=-1
    BB = dist_distr(OA(:,id));
    A = (BB/N)*P;
    Disc_ = h*[1,A]'+b;
    if Disc_<Disc2-epsilon
        Disc2 = Disc_;
        A0 = A;
        id0 = id;
    end
    id = nchoosek_next(id,n1,n);
end
oa = OA(:,id0);
PD = H*[1,A0]'+B;
end
%{
OA = importdata('OA from web/oa.32.9.4.2.a.txt');
n = 5;
[oa1,A1,id1] = MA_from_OA(OA,n);
[oa2,A2,PD2,Disc2,id2] = OA_LowDisc_oalhd(OA,n,'MD');
%}









