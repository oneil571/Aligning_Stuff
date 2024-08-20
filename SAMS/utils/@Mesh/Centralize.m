function [Area,Center] = Centralize(G,scale)
%Centrializes G
%   scale: scale G to a unit 'ScaleArea'

if iscell(G.F)
    error('Not implemented for non-triangular meshes yet');
end
Center = mean(G.V,2);
G.V = G.V-repmat(Center,1,G.nV);

%% Return if no scaling option allowed
if nargin < 2
    return
end

%% Scale mesh if specified
if strcmp(scale,'ScaleArea')
    Area = G.ComputeSurfaceArea;
    G.V = G.V*sqrt(1/Area);
end

end