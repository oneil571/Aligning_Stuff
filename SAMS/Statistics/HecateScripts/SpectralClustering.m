
% CSC - Consistent spectral clustering of surface regions

sqrtInvD(isinf(sqrtInvD)) = 0;
SignVectors = sqrtInvD*U(:, 2:(numEigs+1));
idx = kmeans(SignVectors, numSegments, 'MaxIter', kMeansMaxIter);

%%% TODO: some idx might be +/-Inf, since sqrtInvD might contain +/-Inf
%%% better insert a piece of code here assigning a non-nan label to +/-Inf
%%% points in idx
[InfIdx,~] = find(isinf(SignVectors));
InfIdx = unique(InfIdx);
nVListCumsum = cumsum(diffMatrixSizeList);


for j=1:length(InfIdx)
    IdxJ = find(nVListCumsum>=InfIdx(j),1);
    ValidVList = 1:meshList{IdxJ}.nV;
    IdxOnG = idx(vIdxArray(IdxJ,1):vIdxArray(IdxJ,2));
    ValidVList(IdxOnG == idx(InfIdx(j))) = [];
    tmpDistMatrix = pdist2(meshList{IdxJ}.V(:,InfIdx(j)-vIdxArray(IdxJ,1)+1)'...
        ,meshList{IdxJ}.V(:,ValidVList)');
    [~,minInd] = min(tmpDistMatrix);
    idx(InfIdx(j)) = idx(ValidVList(minInd)+vIdxArray(IdxJ,1)-1);
end
kIdx = idx;
