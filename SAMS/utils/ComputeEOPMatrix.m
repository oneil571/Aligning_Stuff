function Flows = ComputeEOPMatrix(dists,sourceInds,sinkInds)
% Returns the directed flow matrix required for improving pairwise
% correspondence
%
%Input
%dists: matrix of size m x n
%
%Output:
%Flows: cell of size m x n with each entry denoting the directed adjacency
%matrix of flows from i to j

[m,n] = size(dists);
dists = .5*dists+.5*dists';         %symmetrize as sanity check
Flows = cell(length(sourceInds),length(sinkInds));

progressbar
for i = 1:length(sourceInds)
    for j = 1:length(sinkInds)
        fprintf('%d %d \n',i,j);
        Flows{i,j} = sparse(m,n);
        if sourceInds(i) ~= sinkInds(j)
            dummy = sparse(m,n);
            d_i = graphshortestpath(sparse(dists),sourceInds(i));
            d_j = graphshortestpath(sparse(dists),sinkInds(j));
            for k = 1:m
                for q = 1:n
                    if d_i(k) < d_i(q)
                        if d_j(k) > d_j(q)
                            dummy(k,q) = 1;
                        end
                    end
                end
            end
            Flows{i,j} = dummy;
        else
            Flows{i,j}(sourceInds(i),sinkInds(j)) = 1;
        end
    end
    progressbar(i/m);
end

end