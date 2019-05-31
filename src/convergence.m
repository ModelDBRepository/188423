
function c=convergence(filePath)

dirList = dir([filePath 'weight.t=*s.txt']);
for t=1:length(dirList)
    weight = load([filePath dirList(t).name]);
    if t==1
        c = zeros(size(weight,1),length(dirList));
    end
    c(:,t) = mean(abs(weight-(weight>0.5)),2);
end
