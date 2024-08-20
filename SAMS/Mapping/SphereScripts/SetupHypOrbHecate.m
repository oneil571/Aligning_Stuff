%% Redo paths
rmpath(genpath([SAMSPath 'Mapping/external']));
addpath(genpath([SAMSPath 'utils/']));
%% Make directories if needed and load relevant data
orbDataPath = [workingPath 'OrbifoldDataHecate/'];
touch(orbDataPath);
load([workingPath 'Names.mat']);
load([workingPath 'newMeshList.mat']);
load([workingPath 'FinalDists.mat']);

%% Find template, learn features to get a good covering of surface
frechMean = find(min(sum(dists.^2))==sum(dists.^2));
options.numGPLmks = numGPLmksHecate;
options.isDisc = 0;
newMeshList{frechMean}.ComputeAuxiliaryInformation(options);

hecateFeatureInds = [];

switch featureMapHecate
    case 'Gauss'
        hecateFeatureInds = newMeshList{frechMean}.Aux.GaussMaxInds;
    case 'Conf'
        hecateFeatureInds = newMeshList{frechMean}.Aux.ConfMaxInds;
end

%% Select initial features
candidateLmks = newMeshList{frechMean}.Aux.GPLmkInds;
if isempty(hecateFeatureInds)
    hecateFeatureInds = newMeshList{frechMean}.Aux.GPLmkInds(1);
    candidateLmks(1) = [];
elseif length(hecateFeatureInds)>maxInitFeaturesHecate
    hecateFeatureInds = hecateFeatureInds(1:maxInitFeaturesHecate);
end

%% Compute landmarks for matching
startInd = length(hecateFeatureInds)+1;
disp('Computing Landmarks for Hecate Maps');
progressbar
for n = startInd:numLmksHecate
    D = newMeshList{frechMean}.PerformFastMarching(hecateFeatureInds);
    [~,ind] = max(D(candidateLmks));
    hecateFeatureInds = [hecateFeatureInds;candidateLmks(ind)];
    candidateLmks(ind)=[];
    progressbar((n-startInd+1)/(numLmksHecate-startInd+1));
end

%% Write output
disp('Writing meshes for Hecate output');
progressbar
for i = 1:length(Names)
    dirString = orbDataPath;
    if ~exist(dirString,'dir')
        mkdir(dirString);
    end
    newMeshList{i}.Write([dirString Names{i} '.off'],'off',options);
    fid = fopen([dirString Names{i} '.txt'],'w');
    for j = 1:length(hecateFeatureInds)
        fprintf(fid,'%d\n',hecateFeatureInds(j));
    end
    fclose(fid);
    progressbar(i/length(Names));
end
