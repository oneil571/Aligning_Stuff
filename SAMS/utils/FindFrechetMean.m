function fMean = FindFrechetMean(G,distr)
wtA = sparse(pdist2(G.V',G.V').*G.A);
possInds = 1:G.nV;
fVar = Inf;
fMean = 0;
for i = 1:length(possInds)
    [dists,~,~] = graphshortestpath(wtA,possInds(i));
    testVar = dists*distr';
    if testVar < fVar
        fVar = testVar;
        fMean = possInds(i);
    end
end

end