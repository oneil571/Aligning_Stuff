%% Set paths
workingPath = [projectDir 'Default/'];
patchPathPrefix = [workingPath 'Patches/subtaxon/Radius_'];
curKey = 'genus';

%% Get categories from samples
groupsToCompare = Groups(curKey);
unGroups = unique(groupsToCompare);

%Separate unknown samples from known ones, and get list of all known sample
%groups
labelGroups = cell(length(unGroups),1);
numLabels = zeros(length(groupsToCompare),1);
mysteryInds = [];
for l = 1:length(groupsToCompare)
    for ul = 1:length(unGroups)
        if ~isempty(strfind(groupsToCompare{l},unGroups{ul}))
            if strcmp(upper(unGroups{ul}),'Eosimias')
                mysteryInds = [mysteryInds l];
            else
                labelGroups{ul} = [labelGroups{ul} l];
            end
            numLabels(l) = ul;
            break;
        end
    end
end

%% Load meshes and find representatives (Frechet Means)
load([workingPath 'newMeshList.mat']);
load([workingPath 'FinalDists.mat']);

frechetMeshes = cell(size(labelGroups));


clear hFrech


for g = 1:length(frechetMeshes)
    curGroup = labelGroups{g};
    curDists = dists(curGroup,curGroup);
    [~,minInd] = min(sum(curDists.^2));
    frechetMeshes{g} = newMeshList{curGroup(minInd)};
end

%% Mean mesh template + features
aveMeanMeshV = zeros(size(newMeshList{1}.V));
for g = 1:length(frechetMeshes)
aveMeanMeshV = aveMeanMeshV + frechetMeshes{g}.V;
end
aveMeanMeshV = aveMeanMeshV/length(frechetMeshes);
aveMeanMesh = Mesh('VF',aveMeanMeshV,newMeshList{1}.F);
options.isDisc = 1;
SmoothCurvatureFields = 3;          %Amount of smoothing applied to curvatures
ConfMaxLocalWidth = 8;              %Conformal factor maxima radius
GaussMaxLocalWidth = 10;            %Gauss curvature maxima radius
MeanMinLocalWidth = 8;              %Absolute mean curvature minima radius
DNEMaxLocalWidth = 8;               %DNE v1 maxima radius
numGPLmks = 200;
options.ConfMaxLocalWidth = ConfMaxLocalWidth;
options.GaussMaxLocalWidth = GaussMaxLocalWidth;
options.MeanMinLocalWidth = MeanMinLocalWidth;
options.DNEMaxLocalWidth = DNEMaxLocalWidth;
options.SmoothCurvatureFields = SmoothCurvatureFields;
options.numGPLmks = numGPLmks;
options.pointCloud = 0;
[aveMeanMesh,aveMeanMesh.Aux] = aveMeanMesh.ComputeAuxiliaryInformation(options);

%% Get random sample based on conf max
numPseudolmks = 100;
%FPSInds = aveMeanMesh.GeodesicFarthestPointSampling(numPseudolmks,aveMeanMesh.Aux.ConfMaxInds);
FPSInds = unique(aveMeanMesh.GetGPLmk(numPseudolmks));
numPseudolmks = length(FPSInds);
avePtCloud = aveMeanMesh.V(:,FPSInds);
avePtCloud = avePtCloud-mean(avePtCloud,2);
avePtCloud = avePtCloud/norm(avePtCloud,'fro');
% Realign sample
for i = 1:length(newMeshList)
    newMeshList{i}.V = newMeshList{i}.V -mean(newMeshList{i}.V(:,FPSInds),2);
    newMeshList{i}.V = newMeshList{i}.V/norm(newMeshList{i}.V(:,FPSInds),'fro');
    [U,~,V] = svd(newMeshList{i}.V(:,FPSInds)*(avePtCloud'));
    for q = 1:newMeshList{i}.nV
        newMeshList{i}.V(:,q) = V*U'*newMeshList{i}.V(:,q);
    end
end
% Build CSV
load([workingPath 'Names.mat']);

fid = fopen(['MysterySample_' num2str(numPseudolmks) '_GP.csv'],'w');

curLine = [curKey ',SpecimenName,'];
for i = 1:numPseudolmks
    curLine = [curLine 'X ' num2str(i) ',Y ' num2str(i) ',Z ' num2str(i)];
    if i == numPseudolmks
        curLine = [curLine '\n'];
    else
        curLine = [curLine ','];
    end
end
fprintf(fid,curLine);

for n = 1:length(Names)
    curLine = [groupsToCompare{n} ',' Names{n} ','];
    for i = 1:numPseudolmks
        curLine = [curLine num2str(newMeshList{n}.V(1,FPSInds(i)))];
        curLine = [curLine ',' num2str(newMeshList{n}.V(2,FPSInds(i)))];
        curLine = [curLine ',' num2str(newMeshList{n}.V(3,FPSInds(i)))];
        if i == numPseudolmks
            curLine = [curLine '\n'];
        else
            curLine = [curLine ','];
        end
    end
    fprintf(fid,curLine);
end
fclose(fid);   


% And now morphologika
fid = fopen(['MysterySample_' num2str(numPseudolmks) '_Morphologika_GP.txt'],'w');
fprintf(fid,'[Individuals]\n');
fprintf(fid,[num2str(length(Names)) '\n']);
fprintf(fid,'[landmarks]\n');
fprintf(fid,[num2str(numPseudolmks) '\n']);
fprintf(fid,'[dimensions]\n3\n[names]\n');
for n = 1:length(Names)
    fprintf(fid,[Names{n} '\n']);
end
fprintf(fid,'\n[rawpoints]\n\n');
for n = 1:length(Names)
    fprintf(fid,['''' Names{n} '\n\n']);
    for i = 1:numPseudolmks
        curStr = [num2str(newMeshList{n}.V(1,FPSInds(i))) ' ' ...
            num2str(newMeshList{n}.V(2,FPSInds(i))) ' ' ...
            num2str(newMeshList{n}.V(3,FPSInds(i))) '\n'];
        fprintf(fid,curStr);
        if i == numPseudolmks
            fprintf(fid,'\n');
        end
    end
end
fclose(fid);