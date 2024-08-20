function Flows = ComputeDirectedFlows_Row(dists,i)
% Returns the directed flow matrix a la "Eyes on the Prize I"
% ONLY DOES THIS FOR ONE ROW AT A TIME
%Input
%dists: matrix of size m x n
%
%Output:
%Flows: cell of size m x n with each entry denoting the directed adjacency
%matrix of flows from i to j

[m,n] = size(dists);
dists = .5*dists+.5*dists';         %symmetrize as sanity check
Flows = cell(1,n);
    parfor j = 1:n
        fprintf('%d %d \n',i,j);
        Flows{1,j} = sparse(m,n);
        if i ~= j
            for k = 1:m
                for q = 1:n
                    if dists(i,k) < dists(i,q)
                        if dists(j,k) > dists(j,q)
                            Flows{1,j}(k,q) = 1;
                        end
                    end
                end
            end
        else
            Flows{1,j}(i,j) = 1;
        end
    end

end