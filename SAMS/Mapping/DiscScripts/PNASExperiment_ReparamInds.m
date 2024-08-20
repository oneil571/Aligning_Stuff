%% Set definitions
num_iter = 500;
TolX = 1e-10;
TolFun = 1e-6;
load([workingPath 'GPDists.mat']);
load([workingPath 'Names.mat']);
load([workingPath 'MappingData/MatchesPairs_Thresheld.mat']);
addpath(genpath('./Mapping/external/AQP_toolbox/'));
frechMean = find(min(sum(GPDists.^2))==sum(GPDists.^2));

samplesPath = [workingPath 'ProcessedMAT/'];
meshList = cell(length(Names),1);

disp('Loading meshList...')
progressbar
for i = 1:length(Names)
    progressbar(i/length(Names))
    load([samplesPath Names{i} '.mat']);
    meshList{i} = G;
end
meshList = meshList;

TextureCoords1 = cell(length(Names),1);
TextureCoords2 = cell(length(Names),1);
TextureCoordsRev1 = cell(length(Names),1);
TextureCoordsRev2 = cell(length(Names),1);
maps = cell(length(Names),1);
translation = cell(length(Names),1);

%% Compute registrations
disp('Computing registrations');

frechUni = Mesh('VF',meshList{frechMean}.Aux.IsoUniV(1:2,:),meshList{frechMean}.F);
    

progressbar
for i = 1:length(Names)
    disp(['Mesh ' num2str(i) ' out of ' num2str(length(Names))]);
    %progressbar(i/length(Names));
    if i == frechMean
        TextureCoords1{i} = frechUni.V(1:2,:);
        TextureCoords2{i} = TextureCoords1{i};
        TextureCoordsRev1{i} = TextureCoords1{i};
        TextureCoordsRev2{i} = TextureCoords1{i};
        continue;
    end
    curUni = Mesh('VF',meshList{i}.Aux.IsoUniV(1:2,:),meshList{i}.F);
    curPairs = matchesPairs{i};
   

    % global alignment -- apply global rigid tranformation to {1} based on correspondences
    [T.R, T.t] = findBestRigidMotion(curUni.V(1:2,curPairs(:,1))',frechUni.V(1:2,curPairs(:,2))',1);
    % compute linear system (corresponding to landmark matches)
    V = (curUni.V(1:2,:)'+T.t)*T.R';
    if det(T.R)>0
        F =  curUni.F';
    else
        F =  curUni.F([1 3 2],:)';
    end
    [eq_lhs,eq_rhs] = indCoordsToLinearSystem(curUni.V(1:2,:)',...
        curPairs(:,1), ...
        frechUni.V(1:2,curPairs(:,2))');
    % setup optimization problem
    
    optimProblem = OptimProblemIsoDist(V, F, eq_lhs, eq_rhs, [], 5);
    % setup solver
    solverInterp = OptimSolverAcclQuadProx('QP', optimProblem, false, true, true);
    % solve
    logInterp = solverInterp.solveTol(TolX, TolFun ,num_iter);
    % store registered parameterizations
    
    
    %rslt_GP = meshList{i}.InterpolateCorrespondences(meshList{frechMean},options);
    %rslt_cP = meshList{i}.ComputeContinuousProcrustes(meshList{frechMean},options);
    
    %if rslt_GP.cPdist < rslt_cP.cPdist
        TextureCoords1{i} = solverInterp.X';
        TextureCoords2{i} = frechUni.V;
        %maps{i} = rslt_GP.cPmap;
        %translation{i} = rslt_GP.translation;
        %R{i} = rslt_GP.orthogonal;


end
close all


% figure;
% for i = 1:length(Names)
% h(i) = subplot(5,11,i);
% meshList{i}.draw;
% hold on
% end

%% Reparametrize to frechet mean


totalPoints = frechUni.V(1:2,:)';
frechBd = frechUni.FindOrientedBoundaries();
xq = frechUni.V(1,:); yq = frechUni.V(2,:);
pointsToUse = ones(frechUni.nV,1);

for i = 1:length(Names)
    disp(i)
    if i == frechMean
        continue;
    else
        curMesh = Mesh('VF',TextureCoords1{i},meshList{i}.F);
        curBd = curMesh.FindOrientedBoundaries;
        xv = curMesh.V(1,curBd); yv = curMesh.V(2,curBd);
        isInside = inpolygon(xq,yq,xv,yv);
        pointsToUse = pointsToUse .* (isInside');
    end
end
%Extra step: some vertices may only have one neighbor 
%remaining after deletion. Remove those
totalMesh = Mesh('VF',frechUni.V,frechUni.F); 
totalAdj = totalMesh.A;
load([workingPath 'newLmkInds.mat']);
frechLmks = newLmkInds{frechMean};
nullFlag = 1;
while nullFlag
    nullFlag = 0;
    for i = 1:length(pointsToUse)
        if pointsToUse(i)
            curNbr = find(totalAdj(i,:));
            if sum(pointsToUse(curNbr)) <= 1
                pointsToUse(i) = 0;
                nullFlag = 1;
            end
        end
    end
end


lmksToIgnore = find(ismember(frechLmks,find(~pointsToUse)));
disp(['Ignoring ' num2str(length(lmksToIgnore)) ' landmarks']);
disp(lmksToIgnore);
frechLmks(lmksToIgnore) = [];
curFrechLmks = totalMesh.V(:,frechLmks)';
totalMesh.DeleteVertex(find(~pointsToUse));
for i = 1:length(newLmkInds)
    newLmkInds{i}(lmksToIgnore) = [];
end
%curFrechLmks(lmksToIgnore,:) = [];
reparamFrechLmks = knnsearch(totalMesh.V',curFrechLmks);


[~,triAreas] = totalMesh.ComputeSurfaceArea;
totalMesh.Aux.VertArea = (triAreas'*totalMesh.F2V)/3;
%Use new base parametrized meshList to make brand new parametrizations.
%Start with pointLocation for fast computation
newMeshVerts = cell(length(Names),1);
for i = 1:length(Names)
    newMeshVerts{i} = zeros(3,size(totalMesh.V,2));
end
triArray = cell(length(Names),1);
for i = 1:length(Names)
    if i ==frechMean 
        triArray{i} = triangulation(frechUni.F',frechUni.V(1:2,:)');
    else
        triArray{i} = triangulation(meshList{i}.F',TextureCoords1{i}');
    end
end

totalNaNVerts = [];
for i = 1:length(Names)
    [new_ti,new_bc] = triArray{i}.pointLocation(totalMesh.V(1:2,:)');
    
    if i == frechMean
        for j = 1:size(totalMesh.V,2)
            newMeshVerts{i}(:,j) = meshList{i}.V(:,meshList{i}.F(1,new_ti(j))) + ...
            new_bc(j,2)*meshList{i}.V(:,meshList{i}.F(2,new_ti(j))) + ...
            new_bc(j,3)*meshList{i}.V(:,meshList{i}.F(3,new_ti(j)));
           
        end
    end
    for j = 1:size(totalMesh.V,2)
        if sum(isnan(new_bc(j,:))+isnan(new_bc(j,:))+isnan(new_bc(j,:))) == 0
            newMeshVerts{i}(:,j) = new_bc(j,1)*meshList{i}.V(:,meshList{i}.F(1,new_ti(j))) + ...
            new_bc(j,2)*meshList{i}.V(:,meshList{i}.F(2,new_ti(j))) + ...
           new_bc(j,3)*meshList{i}.V(:,meshList{i}.F(3,new_ti(j)));
        end
    end
    new_nan = find(isnan(new_ti));
    totalNaNVerts = [totalNaNVerts;new_nan];
end

nanVerts = unique(totalNaNVerts);

for k = 1:length(nanVerts)
    for i = 1:length(Names)
        BC = triArray{i}.cartesianToBarycentric((1:size(meshList{i}.F,2))',repmat(totalMesh.V(1:2,nanVerts(k))',size(meshList{i}.F,2),1));
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
        newMeshVerts{i}(:,nanVerts(k)) = meshList{i}.V(:,triArray{i}.ConnectivityList(tind,:))*BC';
    end
end

disp('Creating reparametrized meshes...')
unscaledReparamMeshList= cell(length(Names),1);
save([workingPath 'reparamFrechLmks.mat'],'reparamFrechLmks');
for i = 1:length(Names)
    unscaledReparamMeshList{i} = Mesh('VF',newMeshVerts{i},totalMesh.F);
end

save([workingPath 'unscaledReparamMeshList.mat'],'unscaledReparamMeshList');
