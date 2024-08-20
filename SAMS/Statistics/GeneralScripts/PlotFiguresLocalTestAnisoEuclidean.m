close all
adMeanMesh = Mesh('VF',adMean,meanMesh.F);
euMeanMesh = Mesh('VF',euMean,meanMesh.F);
peroMeanMesh = Mesh('VF',peroMean,meanMesh.F);

EuInds = [0,1,2,9,10,11,12,13] + 1;
AdapisInds = [3,4,5,20,21,22,23] + 1;
PerodicticusInds = [6,7,8,14,15,16,17,18,19] + 1;
indList = {AdapisInds,EuInds,PerodicticusInds};

options.mode = 'Native';
genusList = {'Adapis','Eulemur','Perodicticus'};

curPatch=2;
group1 = 2;
group2 = 3;
group1Mean = euMeanMesh;
group2Mean = peroMeanMesh;


%% MDS
Y = mdscale(dists,3,'Criterion','strain');
%% Do rest
curPrefix = 'D://Work/ToothAndClaw/FrechetMean/SpecialPatchesAniso/Patch_';
    
curResults = [curPrefix num2str(curPatch) '_PatchTestResults.mat'];
curMeans = [curPrefix num2str(curPatch) '_Graph.mat'];

load(curResults); load(curMeans);
figure; title('Pairwise Comparison: Frechet Mean');
hold on; curIndFcn = zeros(adMeanMesh.nV,1);
curIndFcn(anisoPatchInds{curPatch}) = 1;

h(1) = subplot(2,2,1);
newMeshList{1}.ViewFunctionOnMesh(curIndFcn,options);
h(2) = subplot(2,2,2);
newMeshList{9}.ViewFunctionOnMesh(curIndFcn,options);

h(3) = subplot(2,2,3);
group1Patch = Mesh('VF',Frechet_Means(:,:,group1),anisoFaces{curPatch});
group2Patch = Mesh('VF',Frechet_Means(:,:,group2),anisoFaces{curPatch});
group1Patch.draw;

displace = group2Patch.V - group1Patch.V;
displaceColor = dot(displace,group1Patch.ComputeNormal)/sqrt(sum(displace.^2));
displaceColor = .5*displaceColor+.5;
displaceColorMat = [];
for k = 1:length(displaceColor)
    displaceColorMat = [displaceColorMat; [1 1 1]*displaceColor(k)];
end
quiver3(group1Patch.V(1,:),group1Patch.V(2,:),group1Patch.V(3,:),...
    displace(1,:),displace(2,:),displace(3,:),'Color','k');

h(4) = subplot(2,2,4);

group2Patch.draw;
quiver3(group2Patch.V(1,:),group2Patch.V(2,:),group2Patch.V(3,:),...
    -displace(1,:),-displace(2,:),-displace(3,:),'Color','k');

format long
disp([genusList{group1} ' vs ' genusList{group2}])
disp(['P-value: ' num2str(pVal(group1,group2))])
Link1 = linkprop(h(1:2), {'CameraUpVector','CameraPosition', 'CameraTarget', 'CameraViewAngle'});
Link2 = linkprop(h(3:4), {'CameraUpVector','CameraPosition', 'CameraTarget', 'CameraViewAngle'});

%% Get patches
patches = cell(length(newMeshList),1);
meanPatch = Mesh('VF',meanMesh.V(:,anisoPatchInds{curPatch}),anisoFaces{curPatch});
meanPatch.V = meanPatch.V - mean(meanPatch.V,2);
meanPatch.V = meanPatch.V/norm(meanPatch.V,'fro');
for i = 1:length(patches)
    patches{i} = Mesh('VF',newMeshList{i}.V(:,anisoPatchInds{curPatch}),anisoFaces{curPatch});
    patches{i}.V = patches{i}.V-mean(patches{i}.V,2);
    patches{i}.V = patches{i}.V/norm(patches{i}.V,'fro');
    [U,~,V] = svd(patches{i}.V*meanPatch.V');
    for q = 1:size(patches{i}.V,2)
        patches{i}.V(:,q) = V*U'*patches{i}.V(:,q);
    end
end

dists = zeros(length(newMeshList),length(newMeshList));
for i = 1:size(dists,1)
    for j = 1:size(dists,2)
        dists(i,j) = norm(patches{i}.V-patches{j}.V,'fro');
    end
end

EucMeanPatch = zeros(size(patches{i}.V));
for i = indList{group2}
    EucMeanPatch = EucMeanPatch + patches{i}.V;
end
EucMeanPatch = EucMeanPatch/length(indList{group2});
EucMeanPatch = Mesh('VF',EucMeanPatch,patches{1}.F);
[U,~,V] = svd(EucMeanPatch.V*group2Patch.V');
for q = 1:size(EucMeanPatch.V,2)
    EucMeanPatch.V(:,q) = V*U'*EucMeanPatch.V(:,q);
end
figure; hold on;
h2(1) = subplot(1,3,1);
group2Patch.draw;
title('Frechet Mean of Patches')
h2(2) = subplot(1,3,2);
EucMeanPatch.draw;
title('Euclidean Mean of Patches')
h2(3) = subplot(1,3,3);
[U,~,V] = svd(meanPatch.V*group2Patch.V');
for q = 1:size(meanPatch.V,2)
    meanPatch.V(:,q) = V*U'*meanPatch.V(:,q);
end
meanPatch.draw;
title('Patch from Global Mean')
Link3 = linkprop(h2, {'CameraUpVector','CameraPosition', 'CameraTarget', 'CameraViewAngle'});
