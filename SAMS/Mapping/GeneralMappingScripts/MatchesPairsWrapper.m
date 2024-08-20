%% Load necessary variables
load([workingPath 'GPDists.mat']);
load([workingPath 'Flags.mat']);
load([workingPath 'Names.mat']);
meshList = cell(length(Names),1);
for i = 1:length(Names)
    load([workingPath 'ProcessedMAT/' Names{i} '.mat']);
    meshList{i}=G;
end
frechMean = find(sum(GPDists.^2)==min(sum(GPDists.^2)));
matchesPairs = cell(length(Names),1);
distList = triu(GPDists); distList = distList(distList>0);
pathWtTemp = 0.25*mean(distList);
startPathWt = exp(-(mean(distList)^2));
%% Run Matches Pairs Procedure
if ~isfield(Flags,'RefinementComputed') || ForceRefinement
    disp('Computing Pairwise Matchings for Mesh:');
    for i = 1:length(Names)
        disp(i)
        if i ~=frechMean
            numLmks = baseLmks;
            curNbrSize = nbrSize;
            while true
                [testMatches,~] = ExtractMatchesPairs(pathWtTemp,startPathWt,pathWtDecr,...
                    i,frechMean,numLmks,curNbrSize,percDecr,minPerc,workingPath);
                testMatches = unique(testMatches,'rows');
                if size(testMatches,1) >= minAlignMatches || numLmks + lmkIter > maxNumLmks
                    if size(testMatches,1) == 0
                        numLmks = baseLmks;
                        curNbrSize = curNbrSize+1;
                    elseif curNbrSize > maxNbrSize
                        error(['Could not detect any correspondences between Frechet mean and \n'...
                            Names{i} '\n Please verify that maxNbrSize is as desired.\n'...
                            'Program will terminate as remaining analysis cannot be trusted.']);
                    else
                        matchesPairs{i} = testMatches;
                        break;
                    end
                else
                    numLmks = numLmks+lmkIter;
                end
            end
        end
    end
    [bd_f,~] = meshList{frechMean}.FindOrientedBoundaries;
    for i = 1:length(Names)
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
    Flags('RefinementComputed') = 1;
    save([workingPath 'Flags.mat'],'Flags');
end
%above save step done in case of errors

%% Threshold matches
load([workingPath '/MappingData/matchesPairs.mat']);
load([workingPath 'MappingData/FeatureMatches.mat']);


frechMesh = meshList{frechMean};
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
        
        featureMatchesPairs{i}(rowsToDel,:) = [];
        oldMatches = [featureMatchesPairs{i};matchesPairs{i}];
        possibleMatches = [matchesPairs{i};featureMatchesPairs{i}];
        newMatches = featureMatchesPairs{i};
        %newMatches = matchesPairs{i}(1,:);
        
        numMatches = size(newMatches,1);
        while numMatches <= maxNumMatches
            if isempty(possibleMatches)
                break;
            end
            totalDists = zeros(size(possibleMatches,1),1);
            [D_cur,~,~] = meshList{i}.PerformFastMarching(newMatches(:,1));
            [D_frech,~,~] = meshList{frechMean}.PerformFastMarching(newMatches(:,2));
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
save([workingPath '/MappingData/MatchesPairs_Thresheld.mat'],'matchesPairs');