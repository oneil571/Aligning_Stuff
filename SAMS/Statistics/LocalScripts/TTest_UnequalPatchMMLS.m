function [f,p]=TTest_UnequalPatch(embedCoords,Group1,Group2,mean1Ind,mean2Ind)
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
for i = 1:length(embedCoords)
    mean1{i} = embedCoords{i}(mean1Ind,:);
    mean2{i} = embedCoords{i}(mean2Ind,:);
    cov1{i} = (embedCoords{i}(Group1,:)-mean1{i})'*(embedCoords{i}(Group1,:)-mean1{i});
    cov1{i} = cov1{i}/(length(Group1)-1);
    cov2{i} = (embedCoords{i}(Group2,:)-mean2{i})'*(embedCoords{i}(Group2,:)-mean2{i});
    cov2{i} = cov2{i}/(length(Group2)-1);
end
%% Computing actual statistic value
t2 = zeros(length(embedCoords),1);
disp('Computing statistic')
for i=1:length(embedCoords)
   t2(i)=(mean2{i}-mean1{i})*((cov1{i}/n1+cov2{i}/n2)\(mean2{i}-mean1{i})');
end

%% Asymptotic approximation
d = length(embedCoords{1}(1,:));
t2 = n1*n2/(n1+n2)*t2;           %T2 statistics, separated only for education
f=(n1+n2-d-1)/(d*(n1+n2-2))*t2;
p=1-fcdf(f, d, n1+n2-d-1 );
