%% Reparametrize to frechet mean
for i = 1:length(Names)
    if i == frechMean
        continue;
    end
    centroid = mean(meshes{i}.V,2);
    meshes{i}.V = meshes{i}.V - repmat(centroid,1,meshes{i}.nV);
    meshes{i}.V = meshes{i}.V/norm(meshes{i}.V,'fro');
    for j = 1:meshes{i}.nV
        meshes{i}.V(:,j) = R{i}'*meshes{i}.V(:,j);
    end
end


totalPoints = meshes{frechMean}.Aux.UniformizationV(1:2,:)';
%totalPoints = uniquetol(totalPoints,1e-3,'ByRows',true);

totalVertList = (1:size(totalPoints,1))';
norm2TotalPoints = totalPoints(:,1).^2 + totalPoints(:,2).^2;
totalToRemove = totalVertList(norm2TotalPoints > maxReparamRadius^2);
totalPoints(totalToRemove,:) = [];

dTri = delaunay(totalPoints);
totalMesh = Mesh('VF',totalPoints',dTri');

vertFaceRing = CORR_compute_vertex_face_ring(totalMesh.F);
vertsToDelete = [];
for i = 1:length(vertFaceRing)
    if isempty(vertFaceRing{i})
        vertsToDelete = [vertsToDelete i];
    end
end

totalMesh.DeleteVertex(vertsToDelete);
%totalMesh.Centralize;
[~,triAreas] = totalMesh.ComputeSurfaceArea;
totalMesh.Aux.VertArea = (triAreas'*totalMesh.F2V)/3;
%Use new base parametrized meshes to make brand new parametrizations.
%Start with pointLocation for fast computation
newMeshVerts = cell(length(Names),1);
for i = 1:length(Names)
    newMeshVerts{i} = zeros(3,size(totalMesh.V,2));
end

triArray = cell(length(Names),1);
for i = 1:length(Names)
    if i ==frechMean 
        triArray{i} = triangulation(meshes{i}.F',meshes{i}.Aux.UniformizationV(1:2,:)');
    else
        triArray{i} = triangulation(meshes{i}.F',TextureCoords1{i}');
    end
end

totalNaNVerts = [];
for i = 1:length(Names)
    if i == 402
        continue;
    end
    [new_ti,new_bc] = triArray{i}.pointLocation(totalMesh.V(1:2,:)');
    for j = 1:size(totalMesh.V,2)
        if sum(isnan(new_bc(j,:))+isnan(new_bc(j,:))+isnan(new_bc(j,:))) == 0
            newMeshVerts{i}(:,j) = new_bc(j,1)*meshes{i}.V(:,meshes{i}.F(1,new_ti(j))) + ...
            new_bc(j,2)*meshes{i}.V(:,meshes{i}.F(2,new_ti(j))) + ...
           new_bc(j,3)*meshes{i}.V(:,meshes{i}.F(3,new_ti(j)));
        end
    end
    new_nan = find(isnan(new_ti));
    totalNaNVerts = [totalNaNVerts;new_nan];
end

nanVerts = unique(totalNaNVerts);

for k = 1:length(nanVerts)
    for i = 1:length(Names)
        BC = triArray{i}.cartesianToBarycentric((1:size(meshes{i}.F,2))',repmat(totalMesh.V(1:2,nanVerts(k))',size(meshes{i}.F,2),1));
        tind = find(all(BC>-3e-1,2));
        if numel(tind)>=1
            smallestBarCoords = min(BC(tind,:)');
            bestInd = find(smallestBarCoords==max(smallestBarCoords));
            bestInd = bestInd(1);
            tind = tind(bestInd);
        else
            disp('Could not find point in triangulation, should never occur');
            rowMins = min(BC,[],2);
            [~,tind] = min(rowMins);
        end
        BC = BC(tind,:);
        newMeshVerts{i}(:,nanVerts(k)) = meshes{i}.V(:,triArray{i}.ConnectivityList(tind,:))*BC';
    end
end

newMeshList= cell(length(Names),1);
for i = 1:length(Names)
    newMeshList{i} = Mesh('VF',newMeshVerts{i},totalMesh.F);
end

%% Renormalize and compute distances based on new alignment
for i = 1:length(newMeshList)
    if i ~= 402
    newMeshList{i}.V = newMeshList{i}.V-mean(newMeshList{i}.V')';
    newMeshList{i}.V = newMeshList{i}.V/norm(newMeshList{i}.V,'fro');
    end
end
save([workingPath 'newMeshList.mat'],'newMeshList');
% Link = linkprop(h, {'CameraUpVector', 'CameraPosition', 'CameraTarget', 'CameraViewAngle'});
% setappdata(gcf, 'StoreTheLink', Link);

dists = zeros(length(Names),length(Names));

%% Computing distances and embeddings
disp('Computing distances and embedding')
for i = 1:length(Names)
    for j = 1:length(Names)
        dists(i,j) = norm(newMeshList{i}.V-newMeshList{j}.V,'fro');
    end
end
[Y,~] = mdscale(dists,3,'Criterion','strain');
save([workingPath 'FinalDists.mat'],'dists'); save([workingPath 'MDSEmbedding.mat'],'Y');

%% Now compute all texture coordinates via composition
disp('Computng all maps between all surfaces via composition')

touch([workingPath 'TextureCoordsSource/']);
touch([workingPath 'TextureCoordsTarget/']);
progressbar

for i = 1:length(Names)
    TextureCoordsSource = cell(length(Names),1);
    TextureCoordsTarget = TextureCoordsSource;
    TextureCoordsSource{i} = TextureCoords1{i};
    TextureCoordsTarget{i} = TextureCoords1{i};
    parfor j =1:length(Names)
        if i ~= j
            [TextureCoordsSource{j},TextureCoordsTarget{j}] ...
                = ComposeTextures([TextureCoords1(i) TextureCoordsRev1(j)]...
                ,[TextureCoords2(i) TextureCoordsRev2(j)]);
        end
        %progressbar(((i-1)*length(Names)+j)/(length(Names)^2));
    end
    progressbar(i/length(Names));
    save([workingPath 'TextureCoordsSource/TextureCoordsSource_' num2str(i) '.mat'],'TextureCoordsSource');
    save([workingPath 'TextureCoordsTarget/TextureCoordsTarget_' num2str(i) '.mat'],'TextureCoordsTarget');
end
