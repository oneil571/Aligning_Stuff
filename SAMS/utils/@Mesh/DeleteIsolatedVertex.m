function dVInds = DeleteIsolatedVertex(G,options)
% This deletes vertices that belong to at most one face, though not the
% faces. This was used primarily for handling topological discs and is
% present in some mapping code but is not actively used.

if nargin<2
    options = [];
end

display = getoptions(options,'display','off');
exclude_boundary = getoptions(options,'exclude_boundary',0);

dVInds = find(logical((sum(G.F2V)==0)+(sum(G.F2V)==1)));

if exclude_boundary==1
    if ~isfield(G,'BV')
        BV = G.FindBoundaries;
    else
        BV = G.BV;
    end
    [~,BVdVInds] = intersect(dVInds,BV);
    dVInds(BVdVInds) = [];
end

if strcmpi(display,'on') && ~isempty(dVInds)
    figure;G.draw();hold on;
    scatter3(G.V(1,dVInds),G.V(2,dVInds),G.V(3,dVInds),30,'g','filled');
end

G.DeleteVertex(dVInds);

end

