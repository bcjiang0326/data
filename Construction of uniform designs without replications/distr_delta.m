function [y,ave_delta] = distr_delta(q,n,prob,flag)
% 20130204 测试 delta/Disc 的分布
% INPUT:
%   q: 各因子水平数
%   n: 设计的行数
%   prob: 概率值，可为vector
%   flag: 0,CD; 1,WD.
% OUTPUT:
%   ID: 输出 prob 分位数，prob 可以为vector

Iter = 100000;
epsilon = 1e-11;

s = length(q);
dd = ones(s,1);
for i=1:s-1
    dd(i)=prod(q(i+1:end));
end

fname = 'delta_distr.txt';
outfile = fopen(fname,'w');
% 首先生成一个随机无重复的矩阵，并生成其附属的变量。生成 z,sigma1,sigma
[D,ID] = rand_U_Type_orth_design(q,n);
if ~flag
    % 选择CD，计算z,sigma1,sigma,disc
    z = zeros(n,s);
    sigma1 = zeros(n,1);
    sigma = zeros(n,n);
    for i = 1:n
        z(i,:) = 0.5*( (D(i,:)+0.5)./q'-0.5 );
    end
    for i = 1:n-1
        for j = i+1:n
            sigma(i,j) = 1/n^2;
            for k = 1:s
                sigma(i,j)=sigma(i,j)*(1+abs(z(i,k))+abs(z(j,k))-abs(z(i,k)-z(j,k)));
            end
            sigma(j,i)=sigma(i,j);
        end 
    end
    for i = 1:n
        sigma1(i) = 2/n;
        for k = 1:s
            sigma1(i) = sigma1(i)*(1+abs(z(i,k))-2*z(i,k)^2);
        end
        sigma(i,i) = prod(1+2*abs(z(i,:)))/n^2;
    end
    % 初始化 Disc (CD)
    Disc=0;
    for i = 1:n-1
        for j = i+1:n
            Disc = Disc + sigma(i,j);
        end
    end
    Disc=(13/12)^s+2*Disc+sum(diag(sigma))-sum(sigma1);
    % 初始化最优解信息
else 
    % 选择WD，计算sigma,disc
    sigma = diag((1.5)^s/n^2*ones(n,1));
    for i = 1:n-1
        for j = i+1:n
            sigma(i,j) = 1/n^2;
            for k = 1:s
                temp = abs(D(i,k)-D(j,k))/q(k);
                sigma(i,j)=sigma(i,j)*(1.5-temp*(1-temp));
            end
            sigma(j,i)=sigma(i,j);
        end 
    end
    % 初始化 Disc (WD)
    Disc=0;
    for i = 1:n-1
        for j = i+1:n
            Disc = Disc + sigma(i,j);
        end
    end
    Disc=-(4/3)^s+2*Disc+sigma(1,1)*n;    
end



%进入迭代  
for Oiter = 1:Iter
    %产生要交换的列及列中两个元素    
    isrepeat = 0;
    k=randi(s);
    i=randi(n);
    j=randi(n);
    while j==i || D(i,k)==D(j,k)    
        j=randi(n);
    end
    if i > j % 强制 i 为较小的那个，后面方便
        temp=i; i=j; j=temp;
    end
    IDi = ID(i) + (D(j,k)-D(i,k))*dd(k); %此处记录下此二数，
    IDj = ID(j) + (D(i,k)-D(j,k))*dd(k); %方便后面更新 ID
    if ~isempty(find(ID==IDi,1)) || ~isempty(find(ID==IDj,1))        
        isrepeat = 1;        
    end        
    %此时首先考虑交换是否会带来重复行
    %若出现重复行，则重新进行交换
    while isrepeat
        isrepeat = 0;
        k=randi(s);
        i=randi(n);
        j=randi(n);
        while j==i || D(i,k)==D(j,k)
            j=randi(n);
        end
        if i > j
            temp=i; i=j; j=temp;
        end
        %检验是否可能出现重复
        IDi = ID(i) + (D(j,k)-D(i,k))*dd(k);
        IDj = ID(j) + (D(i,k)-D(j,k))*dd(k);
        if ~isempty(find(ID==IDi,1)) || ~isempty(find(ID==IDj,1))
            isrepeat = 1;
        end
    end
    
    %计算 delta
    if ~flag
        %对于 CD, 计算 alpha1,alpha2,beta
        alpha1=(1+abs(z(j,k))-2*z(j,k)^2)/(1+abs(z(i,k))-2*z(i,k)^2);
        alpha2=(1+2*abs(z(j,k)))/(1+2*abs(z(i,k)));
        beta = (1+abs(z(j,k))+abs(z(:,k))-abs(z(j,k)-z(:,k)))...
            ./(1+abs(z(i,k))+abs(z(:,k))-abs(z(i,k)-z(:,k)));    
        %计算delta          
        delta = (alpha1-1)*(sigma1(j)/alpha1-sigma1(i))...
            +(alpha2-1)*(sigma(i,i)-sigma(j,j)/alpha2);
        for t = 1:n
            if t ~= i && t ~= j
                delta = delta ...
                    + 2*(beta(t)-1)*(sigma(i,t)-sigma(j,t)/beta(t));
            end
        end
    else
        %对于 CD, 计算 beta
        tempj = abs(D(j,k)-D(:,k))/q(k);
        tempi = abs(D(i,k)-D(:,k))/q(k);
        beta = ( 1.5-tempj.*(1-tempj) )./( 1.5-tempi.*(1-tempi) ); 
        %计算delta          
        delta = 0;
        for t = 1:n
            if t ~= i && t ~= j
                delta = delta ...
                    + (beta(t)-1)*(sigma(i,t)-sigma(j,t)/beta(t));
            end
        end
        delta = delta*2;        
    end
    
    %记录delta/Disc
    if delta > epsilon
        fprintf(outfile,'%.8f %.8f\n',delta/Disc, delta);
    end
    
    % 更新 D,ID,Disc
    temp=D(i,k);D(i,k)=D(j,k);D(j,k)=temp;
    ID(i) = IDi; ID(j) = IDj;
    Disc = Disc + delta;
    
    if ~flag 
        % 若选择 CD, 更新 z,sigma1,diag(sigma),sigma
        temp=z(i,k);z(i,k)=z(j,k);z(j,k)=temp;            
        sigma1(i)=sigma1(i)*alpha1;
        sigma1(j)=sigma1(j)/alpha1;
        sigma(i,i)=sigma(i,i)*alpha2;
        sigma(j,j)=sigma(j,j)/alpha2;
        for t=1:n
            if t~=i && t~=j
                sigma(i,t)=sigma(i,t)*beta(t);
                sigma(j,t)=sigma(j,t)/beta(t);
                sigma(t,i)=sigma(i,t);
                sigma(t,j)=sigma(j,t);
            end
        end
    else
        % 对于 WD, 更新 sigma
        for t=1:n
            if t~=i && t~=j
                sigma(i,t)=sigma(i,t)*beta(t);
                sigma(j,t)=sigma(j,t)/beta(t);
                sigma(t,i)=sigma(i,t);
                sigma(t,j)=sigma(j,t);
            end
        end
    end
end

data = importdata(fname);
y = quantile(data(:,1),prob);
ave_delta = mean(data(:,2));

end
%{
% 测试代码
lv = 3; s1 = 5; s2 = 7; n = 102; prob = [0.05:0.05:0.95]';
[y_CD1,mCD1] = distr_delta(lv*ones(s1,1),n,prob,0);
[y_CD2,mCD2] = distr_delta(lv*ones(s2,1),n,prob,0);
[y_WD1,mWD1] = distr_delta(lv*ones(s1,1),n,prob,1);
[y_WD2,mWD2] = distr_delta(lv*ones(s2,1),n,prob,1);
%}

    