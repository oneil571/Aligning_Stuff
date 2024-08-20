%% Set definitions

load([workingPath 'GPDists.mat']);
load([workingPath 'Names.mat']);
load([workingPath 'MappingData/MatchesPairs_Thresheld.mat']);
frechMean = find(min(sum(GPDists.^2))==sum(GPDists.^2));

samplesPath = [workingPath 'ProcessedMAT/'];
meshes = cell(length(Names),1);

disp('Loading meshes...')
for i = 1:length(Names)
    load([samplesPath Names{i} '.mat']);
    meshes{i} = G;
end
meshList = meshes;
options.FeatureType = 'ConfMax';
options.NumDensityPnts = 100;
options.AngleIncrement = 0.01;
options.NumFeatureMatch = 4;
options.GaussMinMatch = 'off';

TextureCoords1 = cell(length(Names),1);
TextureCoords2 = cell(length(Names),1);
TextureCoordsRev1 = cell(length(Names),1);
TextureCoordsRev2 = cell(length(Names),1);
maps = cell(length(Names),1);
translation = cell(length(Names),1);


%% Compute initial alignments based on landmarks
progressbar
for i = 1:length(Names)
    progressbar(i/length(Names));
    if i == frechMean
        continue;
    end
    if isempty(matchesPairs{i})
        disp('No landmarks matched, continuing');
        continue;
    end
    options.GPMatches = matchesPairs{i};
    curMatchedLmks = matchesPairs{i};
    %Gather points
    ptCloud_1 = meshes{i}.V(:,curMatchedLmks(:,1));
    ptCloud_2 = meshes{frechMean}.V(:,curMatchedLmks(:,2));

    %centralize and normalize
    ptCloud_1 = ptCloud_1 - repmat(mean(ptCloud_1')',1,size(ptCloud_1,2));
    ptCloud_2 = ptCloud_2 - repmat(mean(ptCloud_2')',1,size(ptCloud_2,2));
    ptCloud_1 = ptCloud_1/norm(ptCloud_1,'fro');
    ptCloud_2 = ptCloud_2/norm(ptCloud_2,'fro');

    [U,~,V] = svd(ptCloud_1*(ptCloud_2'));
    for q = 1:meshes{i}.nV
        meshes{i}.V(:,q) = V*U'*meshes{i}.V(:,q);
    end
end
close all
R = cell(length(Names),1);
%% Compute registrations
disp('Computing registrations');
progressbar
for i = 1:length(Names)
    progressbar(i/length(Names));
    if i == frechMean
        TextureCoords1{i} = meshes{i}.Aux.UniformizationV(1:2,:);
        TextureCoords2{i} = TextureCoords1{i};
        TextureCoordsRev1{i} = TextureCoords1{i};
        TextureCoordsRev2{i} = TextureCoords1{i};
        continue;
    end
    if isempty(matchesPairs{i})
        rslt_cP = meshes{i}.ComputeContinuousProcrustesStable_FixedRotation(meshes{frechMean},options);
        TextureCoords1{i} = rslt_cP.TextureCoords1;
        TextureCoords2{i} = rslt_cP.TextureCoords2;
        maps{i} = rslt_cP.cPmap;
        translation{i} = rslt_cP.translation;
        R{i} = rslt_cP.orthogonal;
    end
    options.GPMatches = matchesPairs{i};

    rslt_GP = meshes{i}.InterpolateCorrespondences(meshes{frechMean},options);
    %rslt_cP = meshes{i}.ComputeContinuousProcrustes(meshes{frechMean},options);
    
    %if rslt_GP.cPdist < rslt_cP.cPdist
        TextureCoords1{i} = rslt_GP.TextureCoords1;
        TextureCoords2{i} = rslt_GP.TextureCoords2;
        maps{i} = rslt_GP.cPmap;
        translation{i} = rslt_GP.translation;
        R{i} = rslt_GP.orthogonal;
    %else
    %    TextureCoords1{i} = rslt_cP.TextureCoords1;
    %    TextureCoords2{i} = rslt_cP.TextureCoords2;
    %    maps{i} = rslt_cP.cPmap;
    %    translation{i} = rslt_cP.translation;
    %    R{i} = rslt_cP.orthogonal;
    %end
    
    options.GPMatches = flip(matchesPairs{i},2);

    rslt_GP = meshes{frechMean}.InterpolateCorrespondences(meshes{i},options);
    %rslt_cP = meshes{frechMean}.ComputeContinuousProcrustes(meshes{i},options);
    %if rslt_GP.cPdist < rslt_cP.cPdist
        TextureCoordsRev1{i} = rslt_GP.TextureCoords1;
        TextureCoordsRev2{i} = rslt_GP.TextureCoords2;
    %else
    %    TextureCoordsRev1{i} = rslt_cP.TextureCoords1;
    %    TextureCoordsRev2{i} = rslt_cP.TextureCoords2;
    %end
end
close all


% figure;
% for i = 1:length(Names)
% h(i) = subplot(5,11,i);
% meshes{i}.draw;
% hold on
% end

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
    newMeshList{i}.V = newMeshList{i}.V-mean(newMeshList{i}.V')';
    newMeshList{i}.V = newMeshList{i}.V/norm(newMeshList{i}.V,'fro');
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
