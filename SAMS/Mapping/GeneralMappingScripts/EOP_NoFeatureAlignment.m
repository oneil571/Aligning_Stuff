%% Threshold matches
load([workingPath '/MappingData/matchesPairs.mat']);
load([workingPath 'MappingData/FeatureMatches.mat']);


frechMesh = meshList{frechMean};
frechMeshGraph = graph(sparse(pdist2(frechMesh.V',frechMesh.V').*frechMesh.A));
frechGraphDists = distances(frechMeshGraph);
meanOnFrechGraph = find(min(sum(frechGraphDists.^2)) == sum(frechGraphDists.^2));
[bd_f,~] = frechMesh.FindOrientedBoundaries;
progressbar
for i = 1:length(Names)
    if i ~= frechMean
        curMesh = meshList{i};
        
        [bd_i,~] = meshList{i}.FindOrientedBoundaries;
        rowsToDel = [];
        for j = 1:size(featureMatchesPairs{i})
            if ismember(featureMatchesPairs{i}(j,1),bd_i) || ismember(featureMatchesPairs{i}(j,2),bd_f)
                rowsToDel = [rowsToDel j];
            end
        end
        
        %if size(featureMatchesPairs{i},1) > 0
        %    featureMatchesPairs{i}(rowsToDel,:) = [];
        %    newMatches = featureMatchesPairs{i}(1,:);
        %end
        oldMatches = [matchesPairs{i}];
        possibleMatches = [matchesPairs{i}];
        %possibleMatches(1,:) = [];
        
        %newMatches = matchesPairs{i}(1,:);
        
        numMatches = size(newMatches,1);
        curMeshGraph = graph(sparse(pdist2(meshList{i}.V',meshList{i}.V').*meshList{i}.A));
        curGraphDists = distances(curMeshGraph);
        meanOnCurGraph = find(min(sum(curGraphDists.^2)) == sum(curGraphDists.^2));
        
        D_cur = distances(curMeshGraph,meanOnCurGraph);
        D_frech = distances(frechMeshGraph,meanOnFrechGraph);
        totalDists = zeros(size(possibleMatches,1),1);
        for k = 1:size(possibleMatches,1)
            totalDists(k)= D_cur(possibleMatches(k,1))+D_frech(possibleMatches(k,2));
        end
        [~,nextInd] = min(totalDists);
        nextInd = nextInd(1);   
        newMatches = possibleMatches(nextInd,:);
        possibleMatches(nextInd,:) = [];
        numMatches = 1;
        
        
        while numMatches < maxNumMatches
            if isempty(possibleMatches)
                break;
            end
            totalDists = zeros(size(possibleMatches,1),1);
            
            if numMatches == 1
                D_cur = distances(curMeshGraph,newMatches(:,1));
                D_frech = distances(frechMeshGraph,newMatches(:,2));
            else
                D_cur = min(distances(curMeshGraph,newMatches(:,1)));
                D_frech = min(distances(frechMeshGraph,newMatches(:,2)));
            end
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
    progressbar(i/length(Names));
end
save([workingPath '/MappingData/MatchesPairs_Thresheld_NoAlignment.mat'],'matchesPairs');