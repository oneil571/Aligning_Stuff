vertex_idx_cumsum;
n = length(Names);
nVertList = zeros(n, 1);
vIdxArray = zeros(n+1, 2);
for i=1:n
    nVertList(i) = meshList{i}.nV;
    vIdxArray(i+1, 1) = vIdxArray(i, 2) + 1;
    vIdxArray(i+1, 2) = vIdxArray(i+1, 1) + meshList{i}.nV-1;
end
vIdxArray(1,:) = [];
vIdxCumSum = cumsum(nVertList);