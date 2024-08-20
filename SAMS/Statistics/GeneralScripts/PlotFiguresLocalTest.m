adMeanMesh = Mesh('VF',adMean,meanMesh.F);
euMeanMesh = Mesh('VF',euMean,meanMesh.F);
peroMeanMesh = Mesh('VF',peroMean,meanMesh.F);

genusList = {'Adapis','Eulemur','Perodicticus'};

group1 = 1;
group2 = 2;

curPatch = 5;
curPrefix = 'D://Work/ToothAndClaw/FrechetMean/SpecialPatchesIso/Patch_';
    
curResults = [curPrefix num2str(curPatch) '_PatchTestResults.mat'];
curMeans = [curPrefix num2str(curPatch) '_Graph.mat'];

load(curResults); load(curMeans);


