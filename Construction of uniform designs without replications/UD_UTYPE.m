function [x wd2 info] = UD_UTYPE(Fun,q0,beq,x0,upperbnd,lowerbnd,opt0)
% 2012/4/24 2012/4/25
% ��֦������1.Convex optimization��2.t��3.U-type design
% ˳����ѣ�����ڵݹ�����в���Ҫ���� vec��ȡ����֮���Ǳ��� j
% Input
%       Fun: function_handle, ���� H �� f��0.5*x'*H*x+f'*x, ��ͬ
%            �ľ�����׼���Ӧ��ͬ�� Fun.
%       q0: �������ӵ�ˮƽ�� s-by-1 ����
%       x0: ������һ����ʼ�Ŀ��нⲢ�� bound = 0.5*x'*H*x
%       beq: Լ������ sum(x) = beq ���� beq < m
%       upperbnd: ��ĿǰΪֹ���ҵ�������ֵ�� default = inf
%       lowerbnd: �����½磬δ���ܹ��ﵽ��  default = -inf
%       opt0: 
%           opt0(1):�Ƿ����͹�Ż���֦���� = 1,default 0
%           opt0(2):�Ƿ���� wd2Ls �е� neighborhood ���� ���еײ����㡣�� = 1��default = 0
%           opt0(3):���� t ��֦����ֵ����Χ(0,1]��default = 1 ��ʾ�����á�
% Output
%       x: ������õ����Ž�
%       info: �洢������Ϣ������
%           info(1): �ڵ�һ������ upperbnd �Ĵ���
%           info(2): �ڵڶ������� upperbnd �Ĵ���
%           info(3): �ݹ���ô���
%           info(4): upperbnd ��ʼֵ
%           info(5): upperbnd ��ֵֹ 
%           info(6): upperbnd ����� lowerbnd �����
%                   i.e. (info(5)-lowerbnd)/lowerbnd.
%                   ȡ-1, �� lowerbnd <= 0.    
%           info(7): ���ʵ�Ҷ�ڵ����



narginchk(3, 7);%, nargin, 'struct'));
nargoutchk(0, 3);%, nargout, 'struct'));

% ȫ�ֱ������� ����1
% H: ����������ϼ�Ŀ¼�м���
% q: ��������ˮƽ����� s-by-1 
% TolXInt: �ж�Ϊ��������ֵ����͹�Ż�����
% options: ͹�Ż��еĲ���
% xx: ��ʼ�㣬���洢���ս��
% epsilon: ���㷨�ľ���
% n: �ܵ��������
% m: ����ˮƽ�������
% s: ���Ӹ���
global H q TolXInt options xx epsilon info1 n m s opt;

% ȫ�ֱ������� ����2
% TriesPerLel����¼ÿ�����ӵ�ÿ��ˮƽ���Ѿ�����Ĵ�����sum(q)-by-1 ����
% NoTriesPerLel����¼ÿ�����ӵ�ÿ��ˮƽ�ϰ��ŵ�0�ĸ�����sum(q)-by-1 ����
% URestrict��s-by-4 ����
%       ��һ�У�������ƽ��ÿˮƽ�������������ȡ����i.e.����1�ĸ��� floor(beq./q)
%       �ڶ��У�����������������ƽ����ˮƽ������ mod(beq,q);
%       �����У�������ƽ��ÿˮƽ��ʵ��Ĵ���������ȡ����i.e.����0�ĸ��� floor((m-beq)./q)
%       �����У������Ӳ��������������ƽ����ˮƽ������  mod(m-beq,q);
% Staff: ��ߣ����������� TriesPerLel �зֽ��λ�ã���ʼλ��ǰһ��
% NGetRestr1����¼��ǰ�����������Ѿ�����ƽ����ˮƽ��
% NGetRestr0: ��¼���������в���������Ѿ�����ƽ����ˮƽ��
% ��������� ���� 0 �ĸ������������ ���� 1 �ĸ�����
global TriesPerLel NoTriesPerLel URestrict Staff NGetRestr1 NGetRestr0; % 2012/4/24

INT = 0;

% ��׳�Լ��
m = prod(q0);
if size(q0,2)~=1
    error('Input variable q must be a s-by-1 vector!\n');
end
if any( q0~=round(q0) ) || any( q0 <= zeros(size(q0)) )
    error('The component of input variable must be positive integer!\n');
end
if beq~=round(beq) || beq <= 0 || beq >= m
    error('Input variable beq must be a positive integer and smaller than prod(q0)!\n');
end
if nargin > 3 + INT && ~isempty(x0) && ( sum(x0==1) ~= beq || sum(x0==1) + sum(x0==0)~= m )
    error('The initial input x0 is wrong!\n');
end


% ȫ�ֱ�������ʼֵ ����1
[H,f] = Fun(q0);
q = q0;
TolXInt = 1e-3;
%options = optimset('display','off','LargeScale','on','Algorithm','interior-point-convex','TolCon',1e-2);
options = optimset('display','off','LargeScale','off','Algorithm','interior-point','TolCon',1e-2);
epsilon = 1e-10;
info1 = zeros(7,1);
n = beq;
% m = length(H);m ֮ǰ�Ѹ�ֵ
s = length(q);
if nargin > 6 + INT
    opt = opt0;
    clear('opt0');
end
if nargin <= 6 + INT || isempty(opt) || opt(3) < 0 || opt(3) > 1
    opt = zeros(3,1);
    opt(3) = 1;
end
if nargin > 3 +INT && ~isempty(x0) 
    xx = x0;
    clear('x0');
else 
    xx = (randperm(m))';
    v1 = xx(1:beq);
    v0 = xx(beq+1:end);
    xx(v1) = 1;
    xx(v0) = 0;
    clear('v0','v1');
end 

%*************************
xx(:) = 0;
%*************************

% ���������ͬ�������
if nargin < 5 + INT || isempty(upperbnd)
    upperbnd = inf;
end
info1(4) = upperbnd;
if nargin < 6 + INT || isempty(lowerbnd)
    lowerbnd = -inf;
end
info1(6) = 1; % ��ʼ��Ϊ�� 0, 0 ��ʾ�Ѵﵽ�½�
if lowerbnd <= 0
    info1(6) = -1;
end


% ȫ�ֱ�������ʼֵ ����2
% �˴���ʼ�ĳ��������2012/4/12 �޸���2012/4/24 2012/4/25
TriesPerLel = zeros(sum(q),1);
NoTriesPerLel = zeros(sum(q),1);
Staff = cumsum([0;q(1:s-1)]);
URestrict = zeros(s,2);
URestrict(:,1) = floor(beq./q);
URestrict(:,2) = mod(beq,q);
URestrict(:,3) = floor((m-beq)./q);
URestrict(:,4) = mod(m-beq,q);%����q-URestrict(:,2)
NGetRestr1 = zeros(s,1);
NGetRestr0 = zeros(s,1);
% �˴�֮�ϳ��������2012/4/12 �޸���2012/4/24 2012/4/25




% a recursive function that processes the BB tree 
% �����õݹ鷨���ж������ı�����ʵ�ַ�֦���編�������滮�����)
bgn = 1; % bgn ��ʾ rec_Branchbound() ��ʼѰ�ҿ���ȡ 1 ��λ��
LelCom = ones(s,1) + Staff; % ��¼��ʼѰ�ҿ���ȡ 1 ��λ�� bgn ��Ӧ��ˮƽ���
upperbnd = rec_Branchbound(bgn,LelCom,f,beq,upperbnd,lowerbnd);

% �������
if info1(6) == 0
    % �ﵽ���½�
    info1(5) = lowerbnd;
else    
    info1(5) = upperbnd;
    if info1(6) ~= -1
        info1(6) = (upperbnd-lowerbnd)/lowerbnd;
    end
end
x = xx;
info = info1;
wd2 = -(4/3)^s + 2*info(5)/n^2;
end


%****************************************************************
% a recursive function that processes the BB tree 
% �����õݹ鷨���ж������ı�����ʵ�ַ�֦���編�������滮�����)
%****************************************************************
function upperbnd = rec_Branchbound(bgn,LelCom,f,beq,upperbnd,lowerbnd) 
% �˴����������2012/4/12
global TriesPerLel NoTriesPerLel NGetRestr1 NGetRestr0; 
% URestrict Ҳ��ȫ�ֱ�����JudgeByU1()��JudgeByU0()������
% ����

global H xx epsilon info1 n m s q opt;
global TolXInt options;
info1(3) = info1(3)+1;


if opt(2) %��ʱ�� wd2Ls �е� neighborhood �������뵽��֦������
    % �����ģ�ѽ��� 4 , ���� 1�� 0����Ϊ1 ʱ�����ô˿����   
    if m-bgn+1 <= 4 || min(beq,m-bgn+1-beq) ==1
        v1 = (bgn:bgn+beq-1)'; v0 = (bgn+beq:m)';
        
        
        %**% 2012/5/1 ����˿飬����� m-bgn+1 == 4������ 1 �ĸ���Ϊ 2
        if m-bgn+1 == 4 && beq ==2
            delta = H(v0(1),v0(2)) + sum(f(v0-bgn+1))-...
                H(v1(1),v1(2))-sum(f(v1-bgn+1));
            if delta < 0
                v = v1;
                v1 = v0;
                v0 = v;
                clear('v');
            end
        end        
        %**% 2012/5/1 �����
        
        
        delta = 0; % F(son) - F(father)
        for k = 1:beq % ��ʱ length(v1) == beq
            for l = 1:(m-bgn+1-beq) % ��ʱ length(v0) == m-bgn+1-beq    
                if length(v1) ~= 1    
                    temp_delta = -( sum(H(v1(k),v1)) - 0.5*H(v1(k),v1(k)) ) + ...
                            ( sum(H(v0(l),[v1(1:k-1);v1(k+1:end);v0(l)])) - ...
                            0.5*H(v0(l),v0(l)) );        
                else
                    temp_delta = 0.5*( H(v0(l),v0(l)) - H(v1,v1) );    
                end
                temp_delta = temp_delta - f(v1(k)-bgn+1) + f(v0(l)-bgn+1);
                if temp_delta < delta    
                    delta = temp_delta;    
                    k_ = k;    
                    l_ = l;    
                end    
            end    
        end
        %if son is better then update i, j, v1, v0    
        if delta < -epsilon              
            temp = v1(k_);    
            v1(k_) = v0(l_);    
            v0(l_) = temp;    
        end
        fval = 0.5*sum(sum(H(v1,v1)));
        if ~isempty(f)
            fval = fval + sum(f(v1-bgn+1));
        end
        if fval < upperbnd - epsilon
            % �˴��ǵ�һ�����ܸ��� upperbnd ��λ��
            upperbnd = fval;
            xx(v1) = 1;
            xx(v0) = 0;
            info1(1) = info1(1)+1;
        end
        
        if upperbnd <= epsilon + lowerbnd
            % �ﵽ���½�
            info1(6) = 0;
        end
        info1(7) = info1(7) + 1;
        return;
    end    
elseif beq < 1e-8 % beq = 0,˵���Ѿ��ߵ�Ҷ�ڵ�
    if upperbnd > epsilon
        % �˴�Ҳ�ǵ�һ�����ܸ��� upperbnd ��λ��
        upperbnd = 0;    
        info1(1) = info1(1)+1;    
    end
    if upperbnd <= epsilon + lowerbnd
        % �ﵽ���½�
        info1(6) = 0;
    end
    info1(7) = info1(7) + 1;
    return;
end
            


if opt(1) % ��ʱ��͹�Ż������֦�ֶ�
    vec = (bgn:m);
    % solve the corresponding QIP model with the integarily constraints removed
    [x,fval,exitflag]=quadprog(H(vec,vec),f,[],[],ones(1,length(vec)),beq,zeros(length(vec),1),ones(length(vec),1),[],options);

    % if the solution is not feasible or the value of the objective function is
    % higher than the current upperbnd return with the input intial solution
    if exitflag<=0 || fval >= upperbnd + epsilon
        return;
    end

    % if the integer-constraint variables turned to be integers within the
    % input tolerance , return
    i = find( abs(x-round(x)) > TolXInt , 1 ); 
    if isempty(i)
        if fval < upperbnd - epsilon    % this solution is better than the current solution hence replace 
        % �˴��ǵڶ��ο��ܸ��� upperbnd ��λ��
            xx(vec) = round(x);     
            upperbnd = fval;
            info1(2) = info1(2)+1;
        end
        if upperbnd <= epsilon + lowerbnd
            % �ﵽ���½�
            info1(6) = 0;
        end
        info1(7) = info1(7) + 1;
        return
    end
    clear('x','fval','vec','i','exitflag');
end



% 2012/4/25 �˿������ǰѰ�ҵ�һ���ܹ����� U-type ��Ƶ�λ��
% ��ǰѰ�ҹ����в��ϸ��� NoTriesPerLel �� NGetRestr0
% ���ӿ�������ʱ splt ΪѰ���ĵ�һ���ܹ����� U-type ��Ƶ�λ��
% LelCom Ϊ splt ��Ӧ��ˮƽ���
% splt+1 ��ӦΪ���������⿪ʼѰ�ҵ�һ���ܹ�ȡ1��λ�á�
NGR0 = zeros(s,1); %���� NGetRestr0 �Ļ��˼���
NTPL = zeros(sum(q),1); % NoTriesPerLelCom �Ļ��˼���
for splt = bgn:m
    [Judge1,Rule1] = JudgeByU1(LelCom);
    if Judge1 
        break; % �ҵ�����ȡ 1 ��λ�ã����� 
    end
    % �������˴���˵����ǰ splt ������ȡ 1
    [Judge0,Rule0] = JudgeByU0(LelCom);
    if ~Judge0
        % ���˴���˵����ǰ splt ��Ҳ����ȡ 0, �޽⡣
        % �������λ���֣��ټ�������ȥ�õ��ı�Ȼ���� U-type �ġ�
        % ��ʱӦ�ý����� NGetRestr0 �� NoTriesPerLelCom , �ٷ��ء�
        if splt > bgn
            NoTriesPerLel = NoTriesPerLel-NTPL;
            NGetRestr0 = NGetRestr0 - NGR0;
        end
        return; 
    end
    NoTriesPerLel(LelCom) = NoTriesPerLel(LelCom) + 1;
    NTPL(LelCom) = NTPL(LelCom)+1;
    NGetRestr0(Rule0) = NGetRestr0(Rule0) + 1;
    NGR0(Rule0) = NGR0(Rule0)+1;
    LelCom = updtLc(LelCom);
end
% 2012/4/25�����



% 2012/4/25 �˿����ھ����Ƿ�ִ��ȡ 1 ��֧������
t = min(abs((n-beq+1)/splt-n/m),abs((beq-1)/(m-splt)-n/m));
if t < opt(3)
    C = 0.5*H(splt,splt)+f(splt-bgn+1); % CΪ���� 1 ����������Ҫ��ȥ�ĳ���
    % splt ��ȡ 1 ��֧�������NGetRestr1, TriesPerLel������Ҫ���ˡ�
    NGetRestr1(Rule1) = NGetRestr1(Rule1) + 1;
    TriesPerLel(LelCom) = TriesPerLel(LelCom) + 1;

    % ע��������� splt+1 ��ʼѰ�ң���Ӧ�� LelCom ҲӦΪ splt+1 ��Ӧ�� LelCom.
    upperbnd1 = rec_Branchbound( splt+1, updtLc(LelCom), f(splt-bgn+2:end,1)+H(splt+1:end,splt), beq-1, upperbnd-C, lowerbnd-C );
    if info1(6) == 0
        % ��ȡ���½磬ֱ�Ӹ��� xx ���˳�������Ҫ�ڸ��� upperbnd��
        % ���յ� upperbnd ���� lowerbnd.
        xx(splt) = 1;
        xx(bgn:splt-1) = 0;    
        return;
    end
    if upperbnd1 < upperbnd - C % if the solution was successfull and gives a better upperbnd
        xx(splt) = 1;
        xx(bgn:splt-1) = 0;
        upperbnd = upperbnd1 + C;
    end
    
    
    %���� NGetRestr1, TriesPerLel
    NGetRestr1(Rule1) = NGetRestr1(Rule1) - 1;
    TriesPerLel(LelCom) = TriesPerLel(LelCom) - 1;
end
clear('Rule1');
% 2012/4/25 �����





% 2012/4/25 �˿����ھ����Ƿ�ִ�� 0 ��֧������
% ȡ 0 ��֧֮ǰҪ�ٴν����жϡ��ж� splt ���Ƿ����ȡ 0 ��֧
[Judge0,Rule0] = JudgeByU0(LelCom);
if ~Judge0
    % ���� NGetRestr0, NoTriesPerLel
    if splt > bgn
        NoTriesPerLel = NoTriesPerLel-NTPL;
        NGetRestr0 = NGetRestr0 - NGR0;
    end
    return;
end

t = min(abs((n-beq)/splt-n/m),abs(beq/(m-splt)-n/m));
if t < opt(3) 
    % splt ��ȡ 0 ��֧�������NGetRestr0, NoTriesPerLel������Ҫ���ˡ�
    NoTriesPerLel(LelCom) = NoTriesPerLel(LelCom) + 1;
    NTPL(LelCom) = NTPL(LelCom)+1;
    NGetRestr0(Rule0) = NGetRestr0(Rule0) + 1;
    NGR0(Rule0) = NGR0(Rule0)+1;
    % ע��������� splt+1��ʼѰ�ң���Ӧ�� LelCom ҲӦΪ splt+1 ��Ӧ�� LelCom.
    upperbnd0 = rec_Branchbound( splt+1, updtLc(LelCom), f(splt-bgn+2:end), beq, upperbnd, lowerbnd );
    if info1(6) == 0
        % ��ȡ���½磬ֱ�Ӹ��� xx ���˳�������Ҫ�ڸ��� upperbnd��
        % ���յ� upperbnd ���� lowerbnd.
        xx(bgn:splt) = 0;      
        return;    
    end
    if upperbnd0 < upperbnd % if the solution was successfull and gives a better upperbnd
        xx(bgn:splt) = 0;
        upperbnd = upperbnd0;
    end
    NoTriesPerLel = NoTriesPerLel-NTPL;
    NGetRestr0 = NGetRestr0 - NGR0;
end
    
end




% ***************************************************************
% �˺��������2012/4/12
% �ж��ڵ�ǰ���Ѵ�ȡ 1 �Ƿ�����U����� Judge1=1��ʾ����ȡ1
% Rule1��¼�ڷ��Ѵ���Ӧ����ˮƽ����дﵽ���ȵ���Щˮƽ��Ӧ�����ӣ�
% ������ȡ1����NGetRestr(Rule1) = NGetRestr1(Rule1)+1;
%****************************************************************
function [Judge1,Rule1] = JudgeByU1( LelCom )
global TriesPerLel URestrict NGetRestr1;
if any( TriesPerLel(LelCom) > URestrict(:,1) );
    Judge1 = false;
    Rule1 = [];
    return;
end
Rule1 = TriesPerLel(LelCom) == URestrict(:,1);
% Rule2 = NGetRestr1 >= URestrict(:,2);
% Rule1��Rule2��Ԫ����Ȳ��ҵ���1ʱ�����жϳ� Judge1 = 0
Rule = Rule1 & NGetRestr1 >= URestrict(:,2);
if any( Rule )
    Judge1 = false;
    Rule1 = [];
    return;
end
Judge1 = true;
end

% �ж��ڵ�ǰ���Ѵ�ȡ 0 �Ƿ�����U����� Judge0=1��ʾ����ȡ0
% Rule0��¼�ڷ��Ѵ���Ӧ����ˮƽ����дﵽ���ȵ���Щˮƽ��Ӧ�����ӣ�
% ������ȡ0����NGetRestr0(Rule1) = NGetRestr1(Rule1)+1;
%****************************************************************
function [Judge0,Rule0] = JudgeByU0( LelCom )
global NoTriesPerLel URestrict NGetRestr0;
if any( NoTriesPerLel(LelCom) > URestrict(:,3) );
    Judge0 = false;
    Rule0 = [];
    return;
end
Rule0 = NoTriesPerLel(LelCom) == URestrict(:,3);
% Rule2 = NGetRestr1 >= URestrict(:,4);
% Rule1��Rule2��Ԫ����Ȳ��ҵ���1ʱ�����жϳ� Judge1 = 0
Rule = Rule0 & NGetRestr0 >= URestrict(:,4);
if any( Rule )
    Judge0 = false;
    Rule0 = [];
    return;
end
Judge0 = true;
end


% 2012/4/24 �ó������LelCom����һ��LelCom����update LelCom
function LelCom = updtLc(LelCom)
global q s Staff;
LelCom(s) = LelCom(s)+1;
for t = s-1:-1:1
    if LelCom(t+1)-Staff(t+1)>q(t+1)
        LelCom(t) = LelCom(t)+1;
        LelCom(t+1) = LelCom(t+1)-q(t+1);
    else
        break;
    end
end
end

%{
%2012/4/24 �ó������LelCom��ǰһ��LelCom���� backward LelCom
function LelCom = backLc(LelCom)
global q s Staff;
LelCom(s) = LelCom(s)-1;
for t = s-1:-1:1
    if LelCom(t+1)-Staff(t+1) < 1
        LelCom(t+1) = LelCom(t+1)+q(t+1);
        LelCom(t) = LelCom(t)-1;
    else
        break;
    end
end
end

% 2012/4/24 �ó������jλ�ö�Ӧ������ˮƽ��ϡ�
function LelCom = calcLc(bgn)
global q Staff;
% �˳���������2012/4/12
% LelCom����¼��������λ�ö�Ӧ������ˮƽ���,local,NonPar
% �˿����Ҫ���������Ӧ��ˮƽ���
LelCom = zeros(s,1);
LelCom(s) = mod(bgn,q(s));
if LelCom(s) == 0
    LelCom(s) = q(s);
end
temp = bgn;
for k = s-1:-1:1
    %����ǰt�����ӵ�ˮƽ��ϻ��ֵȼ��࣬temp ��¼�˸õȼ������������
    temp = (temp-LelCom(k+1))/q(k+1)+1;
    LelCom(k) = mod(temp,q(k));
    if LelCom(k)==0
        LelCom(k) = q(k);
    end
end
LelCom = LelCom+Staff;
% �����
end
%}



