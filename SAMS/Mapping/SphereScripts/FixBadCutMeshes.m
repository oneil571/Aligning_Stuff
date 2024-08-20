%% Recompute features...
addpath(genpath([SAMSPath '/utils/']));
rmpath(genpath([SAMSPath 'Matching/external/']));

options.isDisc = Flags('isDisc');
options.ConfMaxLocalWidth = ConfMaxLocalWidth;              
options.GaussMaxLocalWidth = GaussMaxLocalWidth;           
options.MeanMinLocalWidth = MeanMinLocalWidth;              
options.DNEMaxLocalWidth = DNEMaxLocalWidth;               
options.SmoothCurvatureFields = SmoothCurvatureFields;                   
options.numGPLmks = numGPLmks;                    
options.pointCloud = 0;

badMeshInds = [];
output_path = [workingPath 'ProcessedMAT/'];
for i = badMeshList
    G = Mesh('off',[workingPath 'RawOFF/' Names{i} '.off']);
    %progressbar(i/length(Names));
    try
        [G,Aux] =G.ComputeAuxiliaryInformation(options);
        G.Aux = Aux;
        G.Aux.Name = Names{i};
        save([output_path Names{i} '.mat'],'G');
    catch
        badMeshInds = [badMeshInds i];
    end
end

%% The effect of issues should be minor at this point, use prior information.
%Feature matching...

load([workingPath 'MappingData/FeatureMatches.mat']);

load([workingPath 'GPDists.mat']);
frechInd = find(min(sum(GPDists.^2))==sum(GPDists.^2));
load([workingPath 'ProcessedMAT/' Names{frechInd} '.mat']);
frechMean = G;
featureList = cell(length(Names),1);
for i = badMeshList
    load([workingPath 'ProcessedMAT/' Names{i} '.mat']);
    switch featureMap
        case 'Conf'
            featureList{i} =G.Aux.ConfMaxInds;
        case 'Gauss'
            featureList{i} = G.Aux.GaussMaxInds;
        case 'Mean'
            featureList{i} = G.Aux.MeanMinInds;
        case 'DNE'
            featureList{i} = G.Aux.DNEMaxInds;
    end
end

for i = badMeshList
    load([workingPath 'ProcessedMAT/' Names{i} '.mat']);
    disp(['Computing features for ' Names{i}]);
    if i ~= frechMean
        pairMatches = [];
        map_12 = knnsearch(frechMean.V(:,featureList{frechInd})',...
            G.V(:,featureList{i})');
        map_21 = knnsearch(G.V(:,featureList{i})',...
            frechMean.V(:,featureList{frechInd})');
        for j = 1:length(map_12)
            
            %You did this because you wanted to restrict analysis to GP
            %landmarks, so you projected onto GP landmarks because you were
            %concerned that there would be overlap
            if map_21(map_12(j)) == j
                lmk12 = knnsearch(frechMean.V',...
            G.V(:,featureList{i}(j))');
                lmk21 = knnsearch(G.V',...
            frechMean.V(:,featureList{frechInd}(map_12(j)))');
                weightAdj1 = pdist2(G.V',G.V').*G.A;
                weightAdj2 = pdist2(frechMean.V',frechMean.V')...
                    .*frechMean.A;
                [dist12,~,~] = graphshortestpath(weightAdj2,lmk12,...
                    featureList{frechInd}(map_12(j)));
                [dist21,~,~] = graphshortestpath(weightAdj1,lmk21,...
                   featureList{i}(j));
                if dist12 < maxDistTol && dist21 < maxDistTol
                    featureMatchesPairs{i} = [featureMatchesPairs{i};...
                        featureList{i}(j) featureList{frechInd}(map_12(j))];
                end
            end
        end
    end
end
save([workingPath 'MappingData/FeatureMatches.mat'],'featureMatchesPairs');
%% Refine matches
meshList = cell(length(Names),1);
for i = badMeshList
    load([workingPath 'ProcessedMAT/' Names{i} '.mat']);
    meshList{i} = G;
end
distList = triu(GPDists); distList = distList(distList>0);
pathWtTemp = 0.25*mean(distList);
startPathWt = exp(-(mean(distList)^2));

load([workingPath '/MappingData/matchesPairs.mat']);
disp('Computing Pairwise Matchings for Mesh:');
for i = badMeshList
    disp(i)
    if (i ~=frechInd)
        [matchesPairs{i},~] = ExtractMatchesPairsOnFly(pathWtTemp,startPathWt,pathWtDecr,...
            i,frechInd,numGPLmks,nbrSize,percDecr,minPerc,workingPath,maxDepth);
        matchesPairs{i} = unique(matchesPairs{i},'rows');
    end
end
[bd_f,~] = frechMean.FindOrientedBoundaries;
disp('Removing any matches that involve boundaries')
for i = badMeshList
    rowsToDel = [];
    if i ~= frechMean
        [bd_i,~] = meshList{i}.FindOrientedBoundaries;
        for j = 1:size(matchesPairs{i})
            if ismember(matchesPairs{i}(j,1),bd_i) || ismember(matchesPairs{i}(j,2),bd_f)
                rowsToDel = [rowsToDel j];
            end
        end
        matchesPairs{i}(rowsToDel,:) = [];
    end
end
save([workingPath '/MappingData/matchesPairs.mat'],'matchesPairs');
%% Threshold matches
meshList = cell(length(Names),1);
for i = 1:length(Names)
    load([workingPath 'ProcessedMAT/' Names{i} '.mat']);
    meshList{i}=G;
end
frechMesh = frechMean;
for i = badMeshList
    if i ~= frechInd
        curMesh = meshList{i};
        
        [bd_i,~] = meshList{i}.FindOrientedBoundaries;
        rowsToDel = [];
        for j = 1:size(featureMatchesPairs{i})
            if ismember(featureMatchesPairs{i}(j,1),bd_i) || ismember(featureMatchesPairs{i}(j,2),bd_f)
                rowsToDel = [rowsToDel j];
            end
        end
        
        if size(featureMatchesPairs{i},1) > 0
            featureMatchesPairs{i}(rowsToDel,:) = [];
            newMatches = featureMatchesPairs{i}(1,:);
        end
        oldMatches = [featureMatchesPairs{i};matchesPairs{i}];
        possibleMatches = [matchesPairs{i};featureMatchesPairs{i}];
        possibleMatches(1,:) = [];
        
        %newMatches = matchesPairs{i}(1,:);
        
        numMatches = size(newMatches,1);
        while numMatches <= maxNumMatches
            if isempty(possibleMatches)
                break;
            end
            totalDists = zeros(size(possibleMatches,1),1);
            [D_cur,~,~] = meshList{i}.PerformFastMarching(newMatches(:,1));
            [D_frech,~,~] = frechMesh.PerformFastMarching(newMatches(:,2));
            for k = 1:size(possibleMatches,1)
                totalDists(k)= D_cur(possibleMatches(k,1))+D_frech(possibleMatches(k,2));
            end
            [~,nextInd] = max(totalDists);
            nextInd = nextInd(1);
            if min(D_cur(possibleMatches(nextInd,1)),...
                    D_frech(possibleMatches(nextInd,2))) > minMatchDist ...
                    || minMatchDist <= 0
                newMatches = [newMatches;possibleMatches(nextInd,:)];
                numMatches = numMatches+1;
            end
            possibleMatches(nextInd,:) = [];
        end
        matchesPairs{i} = newMatches;
    end 
end
save([workingPath '/MappingData/MatchesPairs_Thresheld.mat'],'matchesPairs');
disp('New matches saved!')