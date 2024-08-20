function [matchedLmks,R]=ExtractMatchesPairsOnFly(temp,wtB,wtDecr,initMesh,...
    finalMesh,landmarks,neighborhoodSize,percDecr,minPerc,workingPath,...
    maxDepth)

% Script for computing the set of matched landmarks based off of forward
% propagation. Works by gradual matching, decreasing amount of sureness
% needed as time goes on in order to capture fuzzier correspondences. Extra
% means that curvature extrema are added to GP landmarks.

% Input:
% temp: bandwidth parameter for kernel.
% wtB: lower bound on minimum weight of path
% wtDecr: weight decrement in case wtB set too high
% initMesh/finalMesh: meshes to match landmarks for, arguments are
% symmetric
% landmarks: number of Gaussian process landmarks used in matching process.
% neighborhoodSize: size of neighborhood of GP landmark to count as match
% minPerc: initial minimum percent of probability mass required for match
% percDecr: decrement parameter for percent matching
% lowBndPerc: final percentage required.


%% Initialize variables, temporarily here until proper form figured out

%% Load necessary variables
    names_path = [workingPath 'Names.mat'];
    dist_path = [workingPath 'GPDists.mat'];
    flows_path = [workingPath 'Flows.mat'];
    samples_path = [workingPath 'ProcessedMAT/'];
%% Don't unnecessarily run
if initMesh == finalMesh
    error('Same mesh, nothing to do')
end

load(names_path);  
load(dist_path);

cPDistances = GPDists;
load(flows_path);

%% Declare and assign global variables to pass to compile results from internal recursion

%Relabel needs to be deprecated

%Number of landmarks to load/consider
global numLandmarks;    
numLandmarks = landmarks;
t =temp;         %temperature parameter for diffusion
       %matrix of cP distances
cPDistances = (cPDistances+cPDistances')/2;
global meshes
meshes = cell(length(Names),1);
%% Load meshes
for i = 1:length(Names)
    load([samples_path Names{i} '.mat']);
    meshes{i} = G;
end

landmarksToTest_1 = meshes{initMesh}.Aux.GPLmkInds(1:landmarks);
landmarksToTest_2 = meshes{finalMesh}.Aux.GPLmkInds(1:landmarks); 
global vertWeight_12;      
vertWeight_12 = zeros(length(landmarksToTest_1),size(meshes{finalMesh}.V,2));
global vertWeight_21;
vertWeight_21 = zeros(length(landmarksToTest_2),size(meshes{initMesh}.V,2));

%% Load weighted flow matrices
curFlow_12 = Flows{initMesh};
curFlow_21 = curFlow_12';
weightedFlows_12 = sparse((cPDistances.^2)).*(curFlow_12);
weightedFlows_21 = sparse((cPDistances.^2)).*(curFlow_21);
for i = 1:size(weightedFlows_12,1)
    for j = 1:size(weightedFlows_12,2)
        if (weightedFlows_12(i,j) ~=0)
            weightedFlows_12(i,j) = exp(-weightedFlows_12(i,j)/t);
           
        end
        if (weightedFlows_21(i,j) ~=0)
            weightedFlows_21(i,j) = exp(-weightedFlows_21(i,j)/t);
        end
    end
end

%% Run diffusion to gather points and distributions
while true
    DepthFirstSearchPlotting_12(weightedFlows_12,initMesh,finalMesh,landmarksToTest_1,1,wtB,1,initMesh,maxDepth);
    DepthFirstSearchPlotting_21(weightedFlows_21,finalMesh,initMesh,landmarksToTest_2,1,wtB,1,finalMesh,maxDepth);
    if norm(vertWeight_12)>0 && norm(vertWeight_21) > 0
        break;
    else %reset and redo computation with smaller weight
        vertWeight_12 = zeros(length(landmarksToTest_1),size(meshes{finalMesh}.V,2));
        vertWeight_21 = zeros(length(landmarksToTest_2),size(meshes{initMesh}.V,2));
        wtB = wtDecr*wtB;
    end
end
disp('Done with initial computations');
%% Make plots
maps_12 = zeros(size(vertWeight_12,2),1);
maps_21 = zeros(size(vertWeight_21,2),1);
for k = 1:numLandmarks
    curColors_12 = vertWeight_12(k,:);
    curColors_21 = vertWeight_21(k,:);
    %nonZer = find(curColors);
    intensity_12 = max(curColors_12);
    intensity_21 = max(curColors_21);
    maps_12 = maps_12 + (curColors_12'/intensity_12);
    maps_21 = maps_21 + (curColors_21'/intensity_21);
    %scatter3(meshes{finalMesh}.V(1,nonZer), ...
    %    meshes{finalMesh}.V(2,nonZer), ...
    %    meshes{finalMesh}.V(3,nonZer),60,[ones(length(nonZer),1) (ones(length(nonZer),1)-(curColors(nonZer)'/intensity)) (ones(length(nonZer),1)-(curColors(nonZer)'/intensity))],'filled');
end
%maps = (maps>0);
%maps = tanh(maps);



%% Find neighborhoods of points
    curMatchedLmks = [];
    while true
        if minPerc < 0.5
            break;
        end
        neighborhoodList_1 = cell(length(landmarksToTest_1),1);
        neighborhoodList_2 = cell(length(landmarksToTest_2),1);
        adjMat_1 = meshes{initMesh}.A; adjMat_2 = meshes{finalMesh}.A;
        neighMat_1 = sparse(meshes{initMesh}.nV,meshes{initMesh}.nV); neighMat_2 = sparse(meshes{finalMesh}.nV,meshes{finalMesh}.nV);
        for q = 0:neighborhoodSize
            neighMat_1 = neighMat_1 + adjMat_1^q;
            neighMat_2 = neighMat_2 + adjMat_2^q;
        end

        for q = 1:length(landmarksToTest_1)
            neighborhoodList_1{q} = find(neighMat_1(landmarksToTest_1(q),:));
        end
        for q = 1:length(landmarksToTest_2)
            neighborhoodList_2{q} = find(neighMat_2(landmarksToTest_2(q),:));
        end
        
        possibleMatches_1 = cell(length(landmarksToTest_1),1);           %possible correspondences for points on first mesh
        possibleMatches_2 = cell(length(landmarksToTest_2),1);             %possible correspondences for points on second mesh
        possibleMatchWeights_1 = cell(length(landmarksToTest_1),1);
        possibleMatchWeights_2 = cell(length(landmarksToTest_2),1);

        for q = 1:length(landmarksToTest_1)
            if ~isempty(curMatchedLmks)
                    if ismember(q,curMatchedLmks(:,1))
                        continue;
                    end
            end
            inds_12 = find(vertWeight_12(q,:));
            for s = 1:length(landmarksToTest_2)
                if ~isempty(curMatchedLmks)
                    if ismember(s,curMatchedLmks(:,2))
                        continue;
                    end
                end
                matchedInds_12 = inds_12(ismember(inds_12,neighborhoodList_2{s}));
                testWeight_12 = sum(vertWeight_12(q,matchedInds_12))/sum(vertWeight_12(q,:));
                if isempty(curMatchedLmks)
                    if (testWeight_12 > minPerc) 
                        possibleMatches_1{q} = [possibleMatches_1{q} s];
                        possibleMatchWeights_1{q} = [possibleMatchWeights_1{q} testWeight_12];
                    end
                else
                    if (testWeight_12 > minPerc) && ~ismember(q,curMatchedLmks(:,1)) && ~ismember(s,curMatchedLmks(:,2))
                        possibleMatches_1{q} = [possibleMatches_1{q} s];
                        possibleMatchWeights_1{q} = [possibleMatchWeights_1{q} testWeight_12];
                    end
                end
            end
        end
        
        for q = 1:length(landmarksToTest_2)
            if ~isempty(curMatchedLmks)
                    if ismember(q,curMatchedLmks(:,2))
                        continue;
                    end
            end
            inds_21 = find(vertWeight_21(q,:));
            for s = 1:length(landmarksToTest_1)
                if ~isempty(curMatchedLmks)
                    if ismember(s,curMatchedLmks(:,1))
                        continue;
                    end
                end
          
                matchedInds_21 = inds_21(ismember(inds_21,neighborhoodList_1{s}));
                testWeight_21 = sum(vertWeight_21(q,matchedInds_21))/sum(vertWeight_21(q,:));
                if isempty(curMatchedLmks)
                    if (testWeight_21) > minPerc 
                        possibleMatches_2{q} = [possibleMatches_2{q} s];
                        possibleMatchWeights_2{q} = [possibleMatchWeights_2{q} testWeight_21];
                    end
                else
                    if (testWeight_21) > minPerc && ~ismember(q,curMatchedLmks(:,2)) && ~ismember(s,curMatchedLmks(:,1))
                        possibleMatches_2{q} = [possibleMatches_2{q} s];
                        possibleMatchWeights_2{q} = [possibleMatchWeights_2{q} testWeight_21];
                    end
                end
            end
        end
    %% Distill correspondences: only consider landmarks if they mutually pair
    %% Create voting list
    possiblePairs = [];
    possibleWeights = [];

    %first side: adds all mutual and possible matches for first landmark
    for q = 1:length(landmarksToTest_1)
        for s = 1:length(possibleMatches_1{q})
            if ismember(q,possibleMatches_2{possibleMatches_1{q}(s)})
                possiblePairs = [possiblePairs; [q possibleMatches_1{q}(s)]];
                matchingInd = find(q==possibleMatches_2{possibleMatches_1{q}(s)});
                newWeight =.5*possibleMatchWeights_1{q}(s)+...
                    .5*possibleMatchWeights_2{possibleMatches_1{q}(s)}(matchingInd);
                possibleWeights = [possibleWeights newWeight];
            end
        end
    end

    
    if ~isempty(possiblePairs)
        [~,orderedIdx] = sort(possibleWeights);
        orderedPairs = possiblePairs(orderedIdx,:);

        %now extract pairs based on stable mutual matching
        trueCorrespondences = [];
        for q = 1:size(orderedPairs,1)
            if q == 1
                trueCorrespondences = orderedPairs(1,:);
            else
                if (~ismember(orderedPairs(q,1),trueCorrespondences(:,1)) && ...
                        ~ismember(orderedPairs(q,2),trueCorrespondences(:,2)))
                    trueCorrespondences = [trueCorrespondences; orderedPairs(q,:)];
                end
            end
        end
        curMatchedLmks = [curMatchedLmks;trueCorrespondences];
    end
    minPerc = minPerc-percDecr;
    end
    if isempty(curMatchedLmks)
        disp('No matches with given parameters found');
        matchedLmks = [];
        R = eye(3);
        return
    end
  
%% Alignment step

%Gather points
ptCloud_1 = meshes{initMesh}.V(:,landmarksToTest_1(curMatchedLmks(:,1)));
ptCloud_2 = meshes{finalMesh}.V(:,landmarksToTest_2(curMatchedLmks(:,2)));

%centralize and normalize
ptCloud_1 = ptCloud_1 - repmat(mean(ptCloud_1')',1,size(ptCloud_1,2));
ptCloud_2 = ptCloud_2 - repmat(mean(ptCloud_2')',1,size(ptCloud_2,2));
ptCloud_1 = ptCloud_1/norm(ptCloud_1,'fro');
ptCloud_2 = ptCloud_2/norm(ptCloud_2,'fro');

[U,~,V] = svd(ptCloud_1*(ptCloud_2'));
R = V*U';
for q = 1:meshes{initMesh}.nV
    meshes{initMesh}.V(:,q) = V*U'*meshes{initMesh}.V(:,q);
end

matchedLmks = [landmarksToTest_1(curMatchedLmks(:,1))' landmarksToTest_2(curMatchedLmks(:,2))'];
fprintf('%d out of %d landmarks in correspondence \n',size(curMatchedLmks,1),landmarks);
empty_1 = 0;
empty_2 = 0;
%% Determine landmarks that can never correspond
minPerc = 0;
possibleMatches_1 = cell(landmarks,1);
possibleMatches_2 = cell(landmarks,2);
for q = 1:landmarks
    inds_12 = find(vertWeight_12(q,:));
    inds_21 = find(vertWeight_21(q,:));
    for s = 1:landmarks
        matchedInds_12 = inds_12(ismember(inds_12,neighborhoodList_2{s}));
        matchedInds_21 = inds_21(ismember(inds_21,neighborhoodList_1{s}));
        if (length(matchedInds_12) > 0)
            possibleMatches_1{q} = [possibleMatches_1{q} s];
        end
        if (length(matchedInds_21)>0)
            possibleMatches_2{q} = [possibleMatches_2{q} s];
        end
    end
end
emptyList_1 = [];
emptyList_2 = [];
for q = 1:landmarks
    if isempty(possibleMatches_1{q})
        for v = 1:landmarks
            if ismember(q,possibleMatches_2{v})
                break;
            elseif v == landmarks
                empty_1=empty_1+1;
                emptyList_1 = [emptyList_1 q];
            end
        end
    end
    if isempty(possibleMatches_2{q})
        for v = 1:landmarks
            if ismember(q,possibleMatches_1{v})
                break;
            elseif v == landmarks
                empty_2=empty_2+1;
                emptyList_2 = [emptyList_2 q];
            end
        end
    end
end

end
function DepthFirstSearchPlotting_12(weightedFlows,init,final,lmks,weight,wtBnd,curDepth,curPath,maxDepth)
%global uniV2 convexComb totalWeight continueCounter
global meshes
global vertWeight_12
global numLandmarks
    nextVerts = find(weightedFlows(init,:));
    
    for i = 1:length(nextVerts)
        if (curDepth == maxDepth) && (nextVerts(i) ~= final)
            continue;
        end
        nextLmks = knnsearch(meshes{nextVerts(i)}.V',meshes{init}.V(:,lmks)');
        nextWeight = weight*weightedFlows(init,nextVerts(i));
        if (nextWeight < wtBnd) 
            continue;
        end
        
        nextPath = [curPath nextVerts(i)];
        if nextVerts(i) == final

            %scatter3(meshes{final}.V(1,nextLmks),meshes{final}.V(2,nextLmks),meshes{final}.V(3,nextLmks),40,repmat([1 1-(nextWeight-wtBnd)/(1-wtBnd) 1-(nextWeight-wtBnd)/(1-wtBnd)],length(nextLmks),1),'filled');
            for k = 1:numLandmarks
                vertWeight_12(k,nextLmks(k)) = vertWeight_12(k,nextLmks(k)) + nextWeight;
            end
        elseif curDepth < maxDepth
            DepthFirstSearchPlotting_12(weightedFlows,nextVerts(i),final,nextLmks,nextWeight,wtBnd,curDepth+1,nextPath,maxDepth);
        end
    end
            
end

function DepthFirstSearchPlotting_21(weightedFlows,init,final,lmks,weight,wtBnd,curDepth,curPath,maxDepth)
%global uniV2 convexComb totalWeight continueCounter
global meshes
global vertWeight_21
global numLandmarks
    nextVerts = find(weightedFlows(init,:));
    for i = 1:length(nextVerts)
        if (curDepth == maxDepth) && (nextVerts(i) ~= final)
            continue;
        end
        
        nextLmks = knnsearch(meshes{nextVerts(i)}.V',meshes{init}.V(:,lmks)');
        nextWeight = weight*weightedFlows(init,nextVerts(i));
        if (nextWeight < wtBnd) 
            continue;
        end
        nextPath = [curPath nextVerts(i)];
        if nextVerts(i) == final
            nextPath
            %scatter3(meshes{final}.V(1,nextLmks),meshes{final}.V(2,nextLmks),meshes{final}.V(3,nextLmks),40,repmat([1 1-(nextWeight-wtBnd)/(1-wtBnd) 1-(nextWeight-wtBnd)/(1-wtBnd)],length(nextLmks),1),'filled');
            for k = 1:numLandmarks
                vertWeight_21(k,nextLmks(k)) = vertWeight_21(k,nextLmks(k)) + nextWeight;
            end
        elseif curDepth < maxDepth
            DepthFirstSearchPlotting_21(weightedFlows,nextVerts(i),final,nextLmks,nextWeight,wtBnd,curDepth+1,nextPath,maxDepth);
        end
    end
            
end