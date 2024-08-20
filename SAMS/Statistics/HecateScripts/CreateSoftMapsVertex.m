%% Load relevant details
load([workingPath 'Names.mat'])
load([workingPath 'newMeshList.mat']);

%% Compute weighted adjacency matrices of everything in a collection
touch([workingPath 'SoftMapsMatrixVertex/'])
wtAdj = cell(length(newMeshList),1);
for i = 1:length(wtAdj)
    wtAdj{i} = pdist2(newMeshList{i}.V',newMeshList{i}.V').*newMeshList{i}.A;
end
for i = 1:length(wtAdj)
    curAdj = wtAdj{i};
    for j = 1:size(curAdj,1)
        curAdj(j,find(curAdj(j,:))) = exp(-(curAdj(j,find(curAdj(j,:))).^2)/fiberEpsVerts);
        curAdj(j,j) = 1;
        curAdj(j,:) = curAdj(j,:)/sum(curAdj(j,:));
    end
    wtAdj{i} = curAdj;
end
        
end
%% Create full soft maps matrix

progressbar
disp('Creating soft maps matrix');
AugKernel12 = cell(length(meshList),1);
AugKernel21 = AugKernel12;
for i = 1:length(meshList)
    G1 = meshList{i};
    softMapsMatrix = cell(length(meshList),1);
    load([workingPath 'TextureCoordsSource/TextureCoordsSource_' num2str(i) '.mat']);
    parfor j = 1:length(meshList)  
        [~,~,AugKernel12{j},~] = MapSoftenKernel(TextureCoordsSource{j}...
            ,meshList{j}.Aux.UniformizationV(1:2,:),...
            meshList{j}.F,G1.V,meshList{j}.V,fiberEps);
        [~,~,AugKernel21{j},~] = MapSoftenKernel(meshList{j}.Aux.UniformizationV(1:2,:)...
            ,TextureCoordsSource{j},G1.F,meshList{j}.V,G1.V,fiberEps);
        softMapsMatrix{j} = max(AugKernel12{j},AugKernel21{j}');
    end
    save([workingPath 'SoftMapsMatrix/SoftMapsMatrix_' num2str(i) '.mat'],'softMapsMatrix');
    progressbar(i/length(meshList));
end


