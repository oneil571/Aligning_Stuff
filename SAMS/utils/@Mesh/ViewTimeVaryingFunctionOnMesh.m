function ViewTimeVaryingFunctionOnMesh(G, color_data,timePause)
%VIEWFUNCTIONONMESH: visualize a function on a triangular mesh
%   This is the full version, the function needs to be defined at all
%   points of the mesh (no interpolation support)
%   The red color stands for positive values, the darker red means the
%   more positive value; the blue color stands for negative value, the
%   darker blue means the more negative value.
%
%   Tingran Gao, Duke University
%   Email: trgao10@math.duke.edu
%   Feb 5, 2015
%   EDITED Robert J Ravier, 7/30/2020
%   Editing done to standardize plots and remove deprecated
%   code. WILL ALWAYS USE 'jet'

%% Quick test to verify that vector is of right dimension
if size(color_data,1) ~= G.nV
    if size(color_data,2) ~= G.nV
        error('Number of function values does not match number of vertices');
    else
        color_data = color_data';
    end
end

%% Renormalize so as to get full color range
color_data = color_data-min(min(color_data));
if max(max(color_data)) > 0
    color_data = color_data/max(max(color_data))*255;
end


figure; hold on;

%% NEED TO FIGURE OUT HOW TO DEAL WITH KEEPING COLORBAR STATIC
for t = 1:size(color_data,2)
    clear gcf
    colormap('jet');
    G.draw(struct('FaceColor', 'interp', 'FaceVertexCData', color_data(:,t),...
        'CDataMapping', 'direct', 'EdgeColor', 'none', 'FaceAlpha', 1,...
        'AmbientStrength',0.3,'SpecularStrength',0.0));
    pause(timePause);
end

hold on;

end
