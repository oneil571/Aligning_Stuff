function [t2,p]=TTest_Standard_PermutationPatch(embedCoords,Group1,Group2,numPerms)
%
% [h p] =hotelT2(v1,v2)
%
% Computes the Hotelling's T-square statistic for given vector measures v1 and v2
% v1,v2    :  vector measures  
%                For instance if v1 is the size of 2562 x 3 x 47, 
%                v1(:,1,12) is the x-coordinate of the 12-th subject.
%                v1 and v2 do not have to have the same number of subjects.
%
%  h : Hotelling's T2 statistic
%  p : The corresponding p-value
%
% (C) Moo K. Chung 2007-2008
%     Department of Biostatistics and Medical Informatics
%     Waisman Laboratory for Brain Imaging
%     University of Wisconsin-Maison
%  
% email://mkchung@wisc.edu

%% Basic setup
n1 = length(Group1); n2 = length(Group2);
totalSamples = [Group1 Group2];

%% Extract means and covariances
disp('Extracting means and covariances')
mean1 = cell(length(embedCoords),1); mean2 = mean1;
cov1 = mean1; cov2 = mean2;
parfor i = 1:length(embedCoords)
    mean1{i} = mean(embedCoords{i}(Group1,:));
    mean2{i} = mean(embedCoords{i}(Group2,:));
    cov1{i} = cov(embedCoords{i}(Group1,:));
    cov2{i} = cov(embedCoords{i}(Group2,:));
end
%% Computing actual statistic value
t2 = zeros(length(embedCoords),1);
parfor i=1:length(embedCoords)
    t2(i)=(mean2{i}-mean1{i})*(((n1-1)*cov1{i}+(n2-1)*cov2{i})/(n1+n2-2)\(mean2{i}-mean1{i})');
end

%% Either sample permutations or generate all if small enough
if nchoosek(n1+n2,n1) <= numPerms
    numPerms = nchoosek((n1+n2),n1);
    permList = nchoosek(1:(n1+n2),n1);
else
    permList = zeros(numPerms,n1);
    for k = 1:numPerms
        permList(k,:) = randperm((n1+n2),n1);
    end  
end
totalInds = [Group1 Group2];
%% Run permutations
permVals = zeros(numPerms,length(embedCoords));
progressbar
for k = 1:numPerms
    V1 = cell(length(embedCoords),1);
    V2 = V1;
    t2 = zeros(1,length(embedCoords));
    curPerm = permList(k,:);
    parfor i=1:length(embedCoords)
        V1{i} = embedCoords{i}(totalInds(curPerm),:);
        V2{i} = embedCoords{i}(totalInds(setdiff(1:(n1+n2),curPerm)),:);
        t2(i)=(mean(V2{i})-mean(V1{i}))*((((n1-1)*cov(V1{i})+(n2-1)*cov(V2{i}))/(n1+n2-2))\(mean(V2{i})-mean(V1{i}))');
    end
    clear V1
    clear V2
    permVals(k,:) = t2;
    progressbar(k/numPerms);
end

%% Compute nonparametric p value
p = zeros(1,length(embedCoords));
for i = 1:length(embedCoords)
    p(i) = (1+length(find(t2(i) < permVals(:,i))))/(1+numPerms);
end
end
