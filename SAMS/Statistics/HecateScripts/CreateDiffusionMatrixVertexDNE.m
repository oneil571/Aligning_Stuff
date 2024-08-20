%% Load relevant details
load([workingPath 'newMeshList.mat']);
load([workingPath 'Names.mat']);

%% Get some notions of epsilons
for i = 1:length(newMeshList)
    if isempty(newMeshList{i})
        continue;
    end
    newMeshList{i}.Centralize('ScaleArea');
    newMeshList{i}.Aux.Name = Names{i};
end
meshList = newMeshList;
load([workingPath 'FinalDists.mat']);
triDists = triu(dists);
baseEps=median(triDists(triDists>0))^2;

disp('Find guess at fiber epsilon')
progressbar
edgeLengths = [];
for i = 1:length(newMeshList)
    curAdj = triu(pdist2(newMeshList{i}.V',newMeshList{i}.V').*newMeshList{i}.A);
    edgeLengths = [edgeLengths reshape(curAdj(curAdj>0),1,length(find(curAdj)))];
    progressbar(i/length(newMeshList));
end
fiberEpsVerts = median(edgeLengths)^2;
%% Compute horizontal nearest neighbors
baseDistMatrix = dists;
H = sparse(length(newMeshList)*newMeshList{1}.nV,length(newMeshList)*newMeshList{1}.nV);
nnList = zeros(BNN,length(newMeshList));
for i = 1:length(newMeshList)
    [~,idx] = sort(baseDistMatrix(:,i)); 
    nnList(:,i) = idx(2:(1+BNN));
end


%% Compute diffusion
nV = newMeshList{1}.nV;
progressbar
disp('Forming Diffusion');
for i = 1:length(newMeshList)
    curAdj = full(pdist2(newMeshList{i}.V',newMeshList{i}.V').*newMeshList{i}.A);
    curAdj(curAdj == 0) = Inf;
    curAdj = exp(-(curAdj.^2)/fiberEpsVerts);
    curAdj = curAdj + eye(newMeshList{1}.nV);
    curAdj = sparse(diag(sum(curAdj,2))\curAdj);
    
    diffInds = find(sum(ismember(nnList,i)));
    for b = diffInds
        H((b-1)*nV+1:b*nV,(i-1)*nV+1:i*nV) = exp(-(dists(i,b)^2)/baseEps)*curAdj;
    end
    progressbar(i/length(newMeshList))
end

%% Symmetrizing for guaranteed diffusion
disp('Symmetrizing...')
progressbar
for i = 1:length(newMeshList)
    for j = i+1:length(newMeshList)
        H((i-1)*nV+1:i*nV,(j-1)*nV+1:j*nV) = ...
            sparse(max(full(H((i-1)*nV+1:i*nV,(j-1)*nV+1:j*nV)),...
            full(H((j-1)*nV+1:j*nV,(i-1)*nV+1:i*nV))));
        H((j-1)*nV+1:j*nV,(i-1)*nV+1:i*nV) = H((i-1)*nV+1:i*nV,(j-1)*nV+1:j*nV);
    end
    progressbar(i/length(newMeshList));
    if rem(i,100) == 0
        save([workingPath 'DiffusionMatrixVertex.mat'],'H');
    end
end

%% Organization to adapt to past methods
GetVertexOrganization;
diffMatrixSize = vIdxCumSum(end);
diffMatrixSizeList = [0; vIdxCumSum];
diffMatrixSizeList(end) = []; % treated as block shifts
%% Saving
disp('Saving Diffusion');
save([workingPath 'DiffusionMatrixVertex.mat'],'H');


