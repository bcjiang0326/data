data = importdata('two-level best results/results_best_0714.txt');
data = data.data;
%for i = 1:size(data,1)
%    fprintf('%d & %d & %.6f & %.6f & %.1f & %.6f & %.6f & %.1f & %.6f & %.6f & %.1f \\\\\n',data(i,:));
%end
niceACD = zeros(size(data,1),1)>1;
niceMCD = zeros(size(data,1),1)>1;
for i = 1:size(data,1)
    if data(i,3)<=data(i,6) && data(i,3)<=data(i,9)
        niceACD(i) = true;
    end
    if data(i,4)<=data(i,7) && data(i,4)<=data(i,10)
        niceMCD(i) = true;
    end
end

data1 = importdata('CD values on web.txt');
data1 = data1.data;
isbetterThanWeb = zeros(size(data,1),1)>1;
for i = 1:size(data,1)
    if data(i,4)<=data1(i,3)
        isbetterThanWeb(i) = true;
    end
end