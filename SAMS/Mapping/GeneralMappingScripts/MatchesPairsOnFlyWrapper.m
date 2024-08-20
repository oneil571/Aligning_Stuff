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
if ~isfield(Flags,'PutativeMatchesComputed') || ForcePutativeMatching
    disp('Computing Pairwise Matchings for Mesh:');
    for i = 1:length(Names)
        disp(i)
        if (i ~=frechMean)
            [matchesPairs{i},~] = ExtractMatchesPairsOnFly(pathWtTemp,startPathWt,pathWtDecr,...
                i,frechMean,numGPLmks,nbrSize,percDecr,minPerc,workingPath,maxDepth);
            matchesPairs{i} = unique(matchesPairs{i},'rows');
        end
    end
    [bd_f,~] = meshList{frechMean}.FindOrientedBoundaries;
    disp('Removing any matches that involve boundaries')
    progressbar
    for i = 1:length(Names)
        progressbar(i/length(Names));
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
    Flags('PutativeMatchesComputed') = 1;
    save([workingPath 'Flags.mat'],'Flags');
end
%above save step done in case of errors

