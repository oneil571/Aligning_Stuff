function frechMean = ComputeFrechetMean(G,distr)
% Compute approximation to Frechet mean for a distribution on a surface
% Not used much
distMat = pdist2(G.V',G.V');
wtA = sparse(distMat.*G.A);

frechVar = inf*ones(G.nV,1);
suppDistr = find(distr);

for i = 1:length(suppDistr)
    [dists,~,~] = graphshortestpath(wtA,suppDistr(i));
    ind = setdiff(1:G.nV,suppDistr);
    dists(ind) = 0;
    frechVar(suppDistr(i)) = sum(dists.*distr);
end

[~,frechMean] = min(frechVar);
end
