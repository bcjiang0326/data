function y = CD2_From_comp(q,id_bar_D)
% 通过 D 的补设计 bar_D 的 CD^2，计算 CD^2(D)。
% n 为 D 的试验次数，则 N-n 为 bar_D 的试验次数。
% Input: 
%       q: m-by-1 vector 各元素为各因子的水平数
%       id_bar_D: bar_D 的 rank vector
% Output:
%       y: = CD^2(D)
[~,D] = compID(id_bar_D,q);
y = CD2_value( D, q );
end


