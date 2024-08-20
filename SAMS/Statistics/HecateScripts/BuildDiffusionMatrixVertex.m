
% DIFFUSION - Build diffusion kernal matrix from distance matrix


%% process base diffusion

load([workingPath 'FinalDists.mat']);
baseDistMatrix = dists;
baseDistMatrix = baseDistMatrix-diag(diag(baseDistMatrix));
n = size(dists,1);

%% Build diffusion maps weight matrix
[sDists,rowNNs] = sort(baseDistMatrix, 2);
sDists = sDists(:,2:(1+BNN));
rowNNs = rowNNs(:,2:(1+BNN));
baseWeights = sparse(repmat((1:n)' ,1 , BNN),rowNNs , sDists, n, n);
baseWeights = min(baseWeights, baseWeights');
for i = 1:n
    sDists(i,:) = baseWeights(i, rowNNs(i,:));
end
sDists = exp(-sDists.^2/baseEps);

%% build augmented diffusion matrix
load([workingPath 'softMapsMatrix.mat']);
diffMatrixSize = vIdxCumSum(end);
diffMatrixSizeList = [0; vIdxCumSum];
diffMatrixSizeList(end) = []; % treated as block shifts
diffMatrixRowIdx = [];
diffMatrixColIdx = [];
diffMatrixVal = [];

cBack = 0;


disp('Constructing diffusion matrix');
progressbar
for j = 1:n
    [~,idx] = sort(dists(:,j)); idx = idx(2:(1+BNN));
    for nns = 1:BNN
        if (sDists(j, nns) == 0)
            continue;
        end
        k = rowNNs(j, nns);
        
        %%% load texture coordinates
 
        AugKernel12 = softMapsMatrix{j, k};

        % Is the next bit meant to be repeated?        
        [rowIdx, colIdx, val] = find(AugKernel12);
        diffMatrixRowIdx = [diffMatrixRowIdx; rowIdx+diffMatrixSizeList(j)];
        diffMatrixColIdx = [diffMatrixColIdx; colIdx+diffMatrixSizeList(k)];
        diffMatrixVal = [diffMatrixVal; sDists(j, nns)*val];

        [rowIdx, colIdx, val] = find(AugKernel12');
        diffMatrixRowIdx = [diffMatrixRowIdx; rowIdx+diffMatrixSizeList(k)];
        diffMatrixColIdx = [diffMatrixColIdx; colIdx+diffMatrixSizeList(j)];
        diffMatrixVal = [diffMatrixVal; sDists(j, nns)*val];
        progressbar(((j-1)*BNN+nns)/(n*BNN));
    end
end

H = sparse(diffMatrixRowIdx,diffMatrixColIdx,diffMatrixVal,diffMatrixSize,diffMatrixSize);

save([workingPath 'DiffusionMatrixVertex.mat'],'H');