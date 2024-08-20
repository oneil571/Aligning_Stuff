%% Load relevant details
load([workingPath 'Names.mat'])
load([workingPath 'newMeshList.mat'])
load([workingPath 'HecateMaps.mat']);
%% Create full soft maps matrix
softMapsMatrix = cell(length(newMeshList),length(newMeshList));
progressbar
disp('Creating soft maps matrix');
AugKernel12 = cell(length(newMeshList),1);
AugKernel21 = AugKernel12;
for i = 1:length(newMeshList)
    G1 = newMeshList{i};
    parfor j = 1:length(newMeshList)  
        [~,~,AugKernel12{j},~] = MapSoftenKernelSphere(map{j,i},...
            newMeshList{j}.F,G1.V,newMeshList{j}.V,fiberEps);
        [~,~,AugKernel21{j},~] = MapSoftenKernelSphere(map{i,j},...
            G1.F,newMeshList{j}.V,G1.V,fiberEps);
        softMapsMatrix{i,j} = max(AugKernel12{j},AugKernel21{j}');
    end
    progressbar(i/length(newMeshList));
end
save([workingPath 'softMapsMatrix.mat'],'softMapsMatrix');


