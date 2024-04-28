%将 2 进制向量转化为起点为 1 的下角标 
function k = idchange(vec)
k = polyval(vec,2)+1; %多项式按照降序排列
end