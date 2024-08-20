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
%% Run Matches Pairs Procedure
pathWtTemp = 0.25*mean(distList);
startPathWt = exp(-2*(mean(distList)^2));

maxDepth = 4;
touch([workingPath 'EOPExperimentData/MappingData/MinPathWt'])

lmkDists = cell(length(Names),1);
progressbar;
for i = 1:length(Names)
    dummyDists = cell(length(Names),1);
    progressbar(i/length(Names))
    parfor j = 1:length(Names)
        disp(j)
        if i ~= j
            dummyDists{j} = sparse(ExtractMatchesPairsDistributions(pathWtTemp,startPathWt,pathWtDecr,...
                i,j,numGPLmks,nbrSize,percDecr,minPerc,workingPath,maxDepth));
        end
    end
    lmkDists{i} = dummyDists;
    save([workingPath 'EOPExperimentData/MappingData/MinPathWt/lmkDists_BD_2p0.mat'],'lmkDists');
end

startPathWt = exp(-0.5*(mean(distList)^2));


lmkDists = cell(length(Names),1);
progressbar;
for i = 1:length(Names)
    dummyDists = cell(length(Names),1);
    progressbar(i/length(Names))
    parfor j = 1:length(Names)
        if i ~= j
            dummyDists{j} = sparse(ExtractMatchesPairsDistributions(pathWtTemp,startPathWt,pathWtDecr,...
                i,j,numGPLmks,nbrSize,percDecr,minPerc,workingPath,maxDepth));
        end
    end
    lmkDists{i} = dummyDists;
    save([workingPath 'EOPExperimentData/MappingData/MinPathWt/lmkDists_BD_0p5.mat'],'lmkDists');
end

%above save step done in case of errors

