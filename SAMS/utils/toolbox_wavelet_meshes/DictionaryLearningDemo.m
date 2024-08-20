disp('Adding path...')
addpath(genpath('/home/grad/rravier/BaseCode/'));
disp('Path added. \n');
load('/home/grad/rravier/BaseCode/samples/PNAS/a10.mat');

options.mode = 'native';
[BV,~] = G.FindOrientedBoundaries;
numPoints = length(BV);
disp('Spreading boundary points');
newNumBound = 180; %%MAKE SURE DIVISIBLE BY 3, needed for regular hexagon
assert(rem(newNumBound,3)==0);

newRadius = .94;
assert(sum(newRadius.^2 > sum(G.Aux.UniformizationV(:,BV).^2)) == 0); %verify radius isn't too big


cornerVerts = newRadius*[cos((0:6)*pi/3);sin((0:6)*pi/3); zeros(1,7)];
[in,on] = inpolygon(G.Aux.UniformizationV(1,:)',G.Aux.UniformizationV(2,:)',...
    cornerVerts(1,:)',cornerVerts(2,:)');
validInds = find(in+on);
interiorPts = G.Aux.UniformizationV(1:2,validInds)';
flatLocations = [cornerVerts(:,1:6) [0;0;0]];
initTri = delaunay(flatLocations(1:2,:)');
curFlatMesh = Mesh('VF',flatLocations,initTri');
a = .001; %init minimum for subdivision
r = 4; %multiplication parameter for subdivision


[Umin,Umax,Cmin,Cmax,Cmean,Cgauss,Normal] = ComputeCurvature(G,options);
totalCurv = Cmin.^2+Cmax.^2;
curFlatVerts = flatLocations(1:2,:)';
curTri = triangulation(initTri,curFlatVerts);

origTri = triangulation(G.F',G.Aux.UniformizationV(1:2,:)'); %needed to create subdivision surfaces

progMeshList = {};
progFlatMeshList = {};
[ti,bc] = origTri.pointLocation(curFlatVerts);
%need to add some step about dealing with nan's, should not happen at this
%step
newMeshVerts = zeros(3,size(curFlatVerts,1));
for j = 1:size(curFlatVerts,1)
    if sum(isnan(bc(j,:))+isnan(bc(j,:))+isnan(bc(j,:))) == 0
        newMeshVerts(:,j) = bc(j,1)*G.V(:,G.F(1,ti(j))) + ...
        bc(j,2)*G.V(:,G.F(2,ti(j))) + bc(j,3)*G.V(:,G.F(3,ti(j)));
    end
end
curMesh = Mesh('VF',newMeshVerts,curTri.ConnectivityList');
progMeshList = [progMeshList curMesh];
progFlatMeshList = [progFlatMeshList curFlatMesh];
flatFaces = initTri;

numMeshesToMake =7;
progFlatMeshList = cell(7,1); progMeshList = cell(7,1);

progFlatMeshList{1} = curFlatMesh;
progMeshList{1} = curMesh;

options.sub_type = 'linear4';
for i=2:numMeshesToMake
    [v,f] = perform_mesh_subdivision(progFlatMeshList{i-1}.V,progFlatMeshList{i-1}.F,1,options);
    progFlatMeshList{i} = Mesh('VF',v,f);
    curFlatVerts = v';
    curFlatTri = f';
    
    [ti,bc] = origTri.pointLocation(curFlatVerts(:,1:2));
    %need to add some step about dealing with nan's, should not happen at this
    %step
    newMeshVerts = zeros(3,size(curFlatVerts,1));
    totalNaNVerts = [];
    for j = 1:size(curFlatVerts,1)
        if sum(isnan(bc(j,:))+isnan(bc(j,:))+isnan(bc(j,:))) == 0
            newMeshVerts(:,j) = bc(j,1)*G.V(:,G.F(1,ti(j))) + ...
            bc(j,2)*G.V(:,G.F(2,ti(j))) + bc(j,3)*G.V(:,G.F(3,ti(j)));
        end
    end
    nanVerts = find(isnan(ti));
    for k = 1:length(nanVerts)
        BC = origTri.cartesianToBarycentric((1:size(G.F,2))',repmat(curFlatVerts(nanVerts(k),:),size(G.F,2),1));
        tind = find(all(BC>-1e-10,2));
        if numel(tind)>=1
            smallestBarCoords = min(BC(tind,:)');
            bestInd = find(smallestBarCoords==max(smallestBarCoords));
            bestInd = bestInd(1);
            tind = tind(bestInd);
        else
            warning('Could not find point in triangulation, should never occur');
            rowMins = min(BC,[],2);
            [~,tind] = min(rowMins);
        end
        BC = BC(tind,:);
        newMeshVerts(:,nanVerts(k)) = G.V(:,origTri.ConnectivityList(tind,:))*BC';
    end
    curMesh = Mesh('VF',newMeshVerts,curFlatTri');
    progMeshList{i} = curMesh;
end

numLmksExtract = 27;

%i=4 as start to guarantee an interesting number of points seeing as we're
%using the butterfly stencil. final chosen as 6 for computational ease
startMesh = 5;
finalMesh = 6;
for i = startMesh:finalMesh
    [progMeshList{i}.Aux.GPLmkInds,~] = progMeshList{i}.GetGPLmk(numLmksExtract);
end

%Essentially proceed by constructing a lot of pieces of the base mesh.
%Counter used to figure out how redundant this will be

boundVertColl = cell(length(progMeshList),1);

%determine what vertices to filter out due to the subdivision stencil.
%Should have 3-neighborhood to deal without boundary bullshit
for i = startMesh:finalMesh
    [BI,~] = progMeshList{i}.FindOrientedBoundaries;
    boundVertColl{i} = BI;
    vring = progMeshList{i}.ComputeVertexRing;
    for j = 1:length(BI)
        boundVertColl{i} = [boundVertColl{i} vring{BI(j)}];
    end
    boundVertColl{i} = unique(boundVertColl{i});
    BI = boundVertColl{i};
    for j = 1:length(BI)
        boundVertColl{i} = [boundVertColl{i} vring{BI(j)}];
    end
    boundVertColl{i} = unique(boundVertColl{i});
    BI = boundVertColl{i};
    for j = 1:length(BI)
        boundVertColl{i} = [boundVertColl{i} vring{BI(j)}];
    end
end
initDict = sparse(progMeshList{startMesh+1}.nV,1);
curDict = initDict;
for i = startMesh:finalMesh-1
    
    pDistMatrix = squareform(pdist(progMeshList{i}.V'));
    weightedAdj = pDistMatrix .* progMeshList{i}.A;
    curLmkInds = progMeshList{i}.Aux.GPLmkInds;
    admissibleLmks = curLmkInds(find(~(ismember(curLmkInds,boundVertColl{i}))));
    %Toss in points
    for j = 1: length(admissibleLmks)
        basisToAdd = zeros(progMeshList{i+1}.nV,1);
        basisToAdd(admissibleLmks)=1;
        curDict = cat(2,curDict,basisToAdd);
    end
    
    %now edges
    for j = 1:length(admissibleLmks)
        for k = (j+1):length(admissibleLmks)
            [~,path,]=graphshortestpath(weightedAdj,admissibleLmks(j),admissibleLmks(k));
            basisToAdd = zeros(progMeshList{i+1}.nV,1);
            basisToAdd(path)=1;
            curDict = cat(2,curDict,basisToAdd);
        end
    end
    
    %and now faces
    disp('Adding Faces...')
    for j = 1:length(admissibleLmks)
        for k = (j+1):length(admissibleLmks)
            for p = (j+2):length(admissibleLmks)
                [~,path12,]=graphshortestpath(weightedAdj,admissibleLmks(j),admissibleLmks(k));
                [~,path23,]=graphshortestpath(weightedAdj,admissibleLmks(k),admissibleLmks(p));
                [~,path31,]=graphshortestpath(weightedAdj,admissibleLmks(p),admissibleLmks(j));
                basisToAdd = zeros(progMeshList{i+1}.nV,1);
                finalPath = [path12 path23(2:end) path31(2:end)];
                [in,on] = inpolygon(progFlatMeshList{i}.V(1,:)',progFlatMeshList{i}.V(2,:)',...
                    progFlatMeshList{i}.V(1,finalPath)',progFlatMeshList{i}.V(2,finalPath)');
                basisToAdd(find(in+on))=1;
                curDict = cat(2,curDict,basisToAdd);
            end
        end
    end 
    %Have everything subdivided now
    if i < finalMesh
        for p = 1:size(curDict,2)
            curFcn = full(curDict(:,p));
            curFcn = perform_wavelet_mesh_transform({progMeshList{i}.V,progMeshList{i+1}.V},...
                {progMeshList{i}.F,progMeshList{i+1}.F},curFcn, -1, options);
            curDict(:,p) = sparse(curFcn);
        end
    end
end

curDict = [curDict sparse(eye(progMeshList{finalMesh}.nV));];
curDict(:,1) = [];

trueDict = sparse(3*progMeshList{finalMesh}.nV,size(curDict,2));

for i = 1:size(trueDict,2)
    curMat = zeros(3,progMeshList{finalMesh}.nV);
    curMat(:,find(curDict(:,i))) = progMeshList{finalMesh}.V(:,find(curDict(:,i)));
    trueDict(:,i) = curMat(:);
end

y = progMeshList{finalMesh}.V(:);
[b,FitInfo] = lasso(trueDict,y);
norm(y-trueDict*b(:,1));
testMesh = Mesh('VF',reshape(trueDict*b(:,1),3,progMeshList{finalMesh}.nV),progMeshList{finalMesh}.F);
close all
figure;
testMesh.draw