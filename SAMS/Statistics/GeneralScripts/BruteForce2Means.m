function [minVal,minGroupA,pVal,HMin] = BruteForce2Means(x,inds1,inds2)

N = size(x,2);
len=length(inds1);
totalComb = combnk(1:N,len);
numPerm = size(totalComb,1);
values = zeros(1,numPerm);
H = zeros(numPerm,1);

%progressbar
parfor i=1:numPerm
    H(i) = max(sum(ismember(totalComb(i,:),inds1))/length(inds1)...
        ,sum(ismember(totalComb(i,:),inds2))/length(inds2));
    values(i) = KMeans_Cost(x,totalComb(i,:));
    %progressbar(i/numPerm);
end
[minVal,minInd] = min(values);
minGroupA = totalComb(minInd,:);
pVal = (1+sum(H>H(minInd)))/(numPerm+1);
HMin = H(minInd);
end