close all

%% Extract weighted average shape (built specifically for talus sample...)
disp('Creating Average Mean');
euMean = zeros(size(newMeshList{1}.V));
adMean = zeros(size(newMeshList{1}.V));
peroMean = zeros(size(newMeshList{1}.V));
EuInds = [0,1,2,9,10,11,12,13] + 1;
AdapisInds = [3,4,5,20,21,22,23] + 1;
PerodicticusInds = [6,7,8,14,15,16,17,18,19] + 1;

for i = EuInds
    euMean = euMean + newMeshList{i}.V;
end
euMean = euMean/length(EuInds);
for i = AdapisInds
    adMean = adMean + newMeshList{i}.V;
end
adMean = adMean/length(AdapisInds);
for i = PerodicticusInds
    peroMean = peroMean + newMeshList{i}.V;
end
peroMean = peroMean/length(PerodicticusInds);
aveMean = (euMean+adMean+peroMean)/3;
meanMesh = Mesh('VF',aveMean,newMeshList{1}.F);

%% Setup landmarks for patches/sampling
disp('Getting landmarks')
isoRadius = 8;
anisoRadius = 5;
numLmks = 10;
lmks = meanMesh.GetGPLmk(numLmks);

isoPatchInds = cell(numLmks,1);
anisoPatchInds = cell(numLmks*(numLmks-1)/2,1);

%% Get radial patches: based on landmark
disp('Drawing Isotropic patches')
totalNumWalksIso = sparse(newMeshList{1}.nV,newMeshList{1}.nV);
for k = 1:isoRadius+1
    totalNumWalksIso = totalNumWalksIso + newMeshList{1}.A^(k1);
end
progressbar
for i = 1:length(isoPatchInds)
    progressbar(i/length(isoPatchInds))
    isoPatchInds{i} = find(totalNumWalksIso(lmks(i),:));
end

%% Get patches based on ridges: ridges defined by geodesics between two lmks
disp('Drawing Anisotropic patches')
weightedGraph = meanMesh.A.*pdist2(meanMesh.V',meanMesh.V');

totalNumWalksAniso = sparse(newMeshList{1}.nV,newMeshList{1}.nV);
progressbar;
for k = 1:anisoRadius+1
    totalNumWalksAniso = totalNumWalksAniso + newMeshList{1}.A^(k-1);
end

cnt = 0;
for i = 1:length(lmks)
    for j = i+1:length(lmks)
        progressbar(cnt/length(anisoPatchInds));
        cnt = cnt + 1;
        [~,curPath,~] = graphshortestpath(weightedGraph,lmks(i),lmks(j));
        curInds = [];
        for k = 1:length(curPath)
            curInds = [curInds find(totalNumWalksAniso(curPath(k),:))];
        end
        anisoPatchInds{cnt} = unique(curInds);
    end
end

%% Setup output
disp('Creating Output Directory')
isoOutputDir = 'D:/Work/ToothAndClaw/FrechetMean/SpecialPatchesIso/';
anisoOutputDir = 'D:/Work/ToothAndClaw/FrechetMean/SpecialPatchesAniso/';
touch(isoOutputDir); touch(anisoOutputDir);

%% Write patches
disp('Writing Isotropic Patches')
progressbar
for i = 1:length(isoPatchInds)
    progressbar(i/length(isoPatchInds));
    disp(i)
    curPatch = zeros(3,length(isoPatchInds{i}),24);
    meanPatch = meanMesh.V(:,isoPatchInds{i});
    meanPatch = meanPatch-mean(meanPatch,2);
    meanPatch = meanPatch/norm(meanPatch,'fro');
    for j = 1:length(newMeshList)
        curPatchJ = newMeshList{j}.V(:,isoPatchInds{i});
        curPatchJ = curPatchJ - mean(curPatchJ,2);
        curPatchJ = curPatchJ/norm(curPatchJ,'fro');
        [U,~,V] = svd(curPatchJ*(meanPatch'));
        for q = 1:size(curPatchJ,2)
            curPatchJ(:,q) = V*U'*curPatchJ(:,q);
        end
        curPatch(:,:,j) = curPatchJ;
    end
    save([isoOutputDir 'Patch_' num2str(i) '.mat'],'curPatch');
end

disp('Writing Anisotropic Patches')
progressbar
for i = 1:length(anisoPatchInds)
    progressbar(i/length(anisoPatchInds));
    curPatch = zeros(3,length(anisoPatchInds{i}),24);
    meanPatch = meanMesh.V(:,anisoPatchInds{i});
    meanPatch = meanPatch-mean(meanPatch,2);
    meanPatch = meanPatch/norm(meanPatch,'fro');
    for j = 1:length(newMeshList)
        curPatchJ = newMeshList{j}.V(:,anisoPatchInds{i});
        curPatchJ = curPatchJ - mean(curPatchJ,2);
        curPatchJ = curPatchJ/norm(curPatchJ,'fro');
        [U,~,V] = svd(curPatchJ*(meanPatch'));
        for q = 1:size(curPatchJ,2)
            curPatchJ(:,q) = V*U'*curPatchJ(:,q);
        end
        curPatch(:,:,j) = curPatchJ;
    end
    save([anisoOutputDir 'Patch_' num2str(i) '.mat'],'curPatch');
end

%% Get Faces
isoFaces = cell(length(isoPatchInds),1); anisoFaces = cell(length(anisoPatchInds),1);

for i = 1:length(isoFaces)
    curCollection = cell(length(newMeshList),1);
    curPatch = isoPatchInds{i};
    faceList = ismember(meanMesh.F,curPatch);
    faceList = find(sum(faceList)>2);
    curFaces = meanMesh.F(:,faceList);
    for j = 1:length(curPatch)
        curFaces(curFaces == curPatch(j)) = j;
    end
    isoFaces{i} = curFaces;
end
for i = 1:length(anisoFaces)
    curCollection = cell(length(newMeshList),1);
    curPatch = anisoPatchInds{i};
    faceList = ismember(meanMesh.F,curPatch);
    faceList = find(sum(faceList)>2);
    curFaces = meanMesh.F(:,faceList);
    for j = 1:length(curPatch)
        curFaces(curFaces == curPatch(j)) = j;
    end
    anisoFaces{i} = curFaces;
end