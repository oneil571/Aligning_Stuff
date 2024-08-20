adMeanMesh = Mesh('VF',adMean,meanMesh.F);
euMeanMesh = Mesh('VF',euMean,meanMesh.F);
peroMeanMesh = Mesh('VF',peroMean,meanMesh.F);

options.mode = 'Native';
genusList = {'Adapis','Eulemur','Perodicticus'};

group1 = 1;
group2 = 3;
group1Mean = adMeanMesh;
group2Mean = peroMeanMesh;

curPatch = 2;
curPrefix = 'D://Work/ToothAndClaw/FrechetMean/SpecialPatchesAniso/Patch_';
    
curResults = [curPrefix num2str(curPatch) '_PatchTestResults.mat'];
curMeans = [curPrefix num2str(curPatch) '_Graph.mat'];

load(curResults); load(curMeans);
figure; title('Pairwise Comparison: Frechet Mean');
hold on; curIndFcn = zeros(adMeanMesh.nV,1);
curIndFcn(anisoPatchInds{curPatch}) = 1;

h(1) = subplot(2,2,1);
group1Mean.ViewFunctionOnMesh(curIndFcn,options);
h(2) = subplot(2,2,2);
group2Mean.ViewFunctionOnMesh(curIndFcn,options);

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
