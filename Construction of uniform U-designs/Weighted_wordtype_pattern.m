function [WB,aveDisc,A,BB,VEC] = Weighted_wordtype_pattern(D,q,DiscMeasure)
%20150621 by Xiaopang
%�ο� Jiang and Ai (2015) Construction of uniform U designs
%1.���� D ���ض� discrepancy �յ��� weighted wordtype pattern
%2.���� D ���ɵ� U design �� ƽ��������
%      ��� D �ĸ����Ӳ��ǰ���ˮƽ���������У����������� D �ĸ��У�
%      ʹ�ø���ˮƽ��Ϊ��������
%INPUT:
%   D: һ���ǶԳ� OA
%   q: n-by-1 ��������¼OA���е�ˮƽ��
%   DiscMeasure: 'CD','WD','MD'
%OUTPUT:
%   WB: discrepancy induced weighted wordtype pattern
%   aveDisc: ƽ��������
%   A: generalized wordlength pattern
%   BB: generalized wordtype pattern
%   VEC: ����������������ÿ��Ԫ�صĸ�ˮƽ���ӵĸ���
[N,n] = size(D);
if nargin < 3
    DiscMeasure = 'CD';
end
if nargin < 2 || isempty(q)
    q = zeros(n,1);
    for i = 1:n
        q(i) = length(unique(D(:,i)));
    end
end
if size(q,2)~=1 || size(q,1)~=n
    error('q must be an n-by-1 vector!\n');
end
if ~issorted(q)
    [q,id] = sort(q);
    D = D(:,id);
end
if ~strcmp(DiscMeasure,'CD') && ~strcmp(DiscMeasure,'WD') && ~strcmp(DiscMeasure,'MD')
    error('DiscMeasure must be CD,WD or MD!\n');
end

uq = unique(q); %���г��ֵĲ�ͬˮƽ��
m = histc(q,[uq-0.5;inf]); 
m(end) = []; %��ˮƽ���ӵĸ�����m(i) Ϊˮƽ��Ϊ uq(i) �����Ӹ���

if strcmp(DiscMeasure,'CD')
    if mod(N,2)==0
        if nargout > 1
            alpha = (13/12)^n-2*(13/12+1/(24*N^2))^n+1.25^n/N...
                -prod((1.25-1./(6*uq)-1/(6*N)).^m)/N;
            a1 = 13/12-1./(6*N*uq);
        end
        a2 = (2*N-2)./(13*N*uq-2);
    else
        if nargout > 1
            alpha = (13/12)^n-2*(13/12-1/(12*N^2))^n+(1.25-1/(4*N^2))^n/N...
                -prod((1.25-1./(6*uq)-1/(6*N)-1/(4*N^2)).^m)/N;
            a1 = 13/12-1./(6*N*uq)-1/(4*N^2);
        end
        a2 = (2*N-2)./(13*N*uq-2-3*uq/N);
    end            
elseif strcmp(DiscMeasure,'WD')
    if nargout > 1
        alpha = -(4/3)^n+1.5^n/N-prod((1.5-1./(3*uq)-1/(3*N)+1./(6*uq.^2)+1./(6*N*uq)).^m)/N;
        a1 = 4/3-1./(3*N*uq)-1/(6*N^2)+1./(6*N*uq.^2)+1./(6*N^2*uq);
    end
    a2 = (N^2*uq-N^2-2*N*uq+N+uq)./(8*N^2*uq.^2-2*N*uq-uq.^2+N+uq);
elseif strcmp(DiscMeasure,'MD')
    if mod(N,2)==0
        if nargout > 1
            alpha = (19/12)^n-2*(19/12+1/(48*N^2))^n+1.75^n/N...
                -prod((1.75-1./(4*uq)-1/(4*N)+1./(12*uq.^2)+1./(12*N*uq)).^m)/N;
            a1 = 19/12-1./(4*N*uq)+1./(12*N*uq.^2)+1./(12*N^2*uq)-1/(12*N^2);
        end
        a2 = (4*N^2*uq-2*N^2-6*N*uq+2*N+2*uq.^2)./(38*N^2*uq.^2-6*N*uq+2*N+2*uq-2*uq.^2);
    else
        if nargout > 1
            alpha = (19/12)^n-2*(19/12+1/(12*N^2))^n+(1.75+1/(8*N^2))^n/N...
                -prod((1.75-1./(4*uq)-1/(4*N)+1./(12*uq.^2)+1./(12*N*uq)+1/(8*N^2)).^m)/N;
            a1 = 19/12-1./(4*N*uq)+1./(12*N*uq.^2)+1./(12*N^2*uq)+1/(24*N^2);
        end
        a2 = (4*N^2*uq-2*N^2-6*N*uq+2*N+2*uq.^2)./(38*N^2*uq.^2-6*N*uq+2*N+2*uq+uq.^2);
    end    
end
BB = genWordTypePtn(D,q);

coef = zeros(1,length(BB));
typeLen = zeros(1,length(BB));
if nargout > 4
    VEC = Id2Design(m+1,(0:length(BB)-1)');
end
for k = 1:length(BB)
    vec = Id2Design(m+1,k-1);
    coef(k) = prod(a2'.^vec);
    typeLen(k) = sum(vec);
end
WB = zeros(1,n);
if nargout > 2
    A = zeros(1,n);
end
for j = 1:n
    id = typeLen==j;
    mj = sum(id);
    Cj = sum(coef(id));
    WB(j) = coef(id)*BB(id)*mj/Cj;
    if nargout > 2
        A(j) = sum(BB(id));
    end
end
if nargout>1
    aveDisc = alpha+prod(a1.^m)*coef*BB;
end
end
%{
oa = importdata('OA from web/MA.36.3.12.2.11.txt');
oa = oa(:,10:14);
[N,n] = size(oa);
q = zeros(n,1);
for i = 1:n
    q(i) = length(unique(oa(:,i)));
end
[WB,y,A,BB,VEC] = Weighted_wordtype_pattern(oa,q,'CD');
%}
