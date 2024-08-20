function [V,F] = PerformGeodesicDelauneyTriangulation(G, VertSampInd, options)
%This will attempt to retriangulate a mesh but has not been used in
%pipeline for a while.

options.method = getoptions(options, 'method', 'slow');
options.verb = getoptions(options, 'verb', 0);

[~,~,Q] = G.PerformFastMarching(VertSampInd);
n = length(G.V);
m = length(VertSampInd);
VQ = Q(G.F); VQ = sort(VQ,1);
VQ = unique(VQ', 'rows')';
d = 1 + (VQ(1,:)~=VQ(2,:)) + (VQ(2,:)~=VQ(3,:));
I = find(d==3); I = sort(I);
z = zeros(n,1);
z(VertSampInd) = (1:m)';
F = z(VQ(:,I));
V = G.V(:,VertSampInd);
F = G.PerformFacesReorientation(V, F, options);


end

