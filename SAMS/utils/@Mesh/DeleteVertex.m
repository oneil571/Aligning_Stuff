function DeleteVertex(G,dVInds)
%Deletes vertices

nVOrigInds = 1:G.nV;
nVOrigInds(dVInds) = Inf;
nVNewInds = cumsum(1-isinf(nVOrigInds));
[dFInds,~] = find(G.F2V(:,dVInds));

G.V(:,dVInds) = [];
G.F(:,dFInds) = [];
G.F = nVNewInds(G.F);

G.F2V = G.ComputeF2V;
G.V2E = G.ComputeV2E;
[G.A,G.E] = G.ComputeAdjacencyMatrix;
G.nV = size(G.V,2);
G.nF = size(G.F,2);
G.nE = size(G.V2E,2);

end

