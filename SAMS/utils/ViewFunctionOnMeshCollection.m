function [h,cbh,Link]=ViewFunctionOnMeshCollection(meshList, fcnList,numPerRow,colorbarBool)
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
if nargin < 3
    error('Not enough arguments')
elseif nargin < 4
    colorbarBool = 1;
end
if length(meshList) ~= length(fcnList)
    error('Number of meshes and functions not equal, please fix');
end
totalFcn = [];
for i = 1:length(fcnList)
    if size(fcnList{i},1) ~= meshList{i}.nV
        if size(fcnList{i},2) ~= meshList{i}.nV
            error('Number of function values does not match number of vertices');
        else
            fcnList{i} = fcnList{i}';
        end
    end
    totalFcn = [totalFcn;fcnList{i}];
end

for i = 1:length(fcnList)
    fcnList{i} = (fcnList{i}-min(totalFcn));
    if max(totalFcn) > 0
        fcnList{i} = fcnList{i}/max(totalFcn)*255;
    end
end


colormap('jet');

numRows = ceil(length(meshList)/numPerRow);


for i = 1:length(meshList)
    h(i) = subplot(numRows,numPerRow,i);
    meshList{i}.draw(struct('FaceColor', 'interp', 'FaceVertexCData', fcnList{i},...
        'CDataMapping', 'direct', 'EdgeColor', 'none', 'FaceAlpha', 1,...
        'AmbientStrength',0.3,'SpecularStrength',0.0));
    if colorbarBool
        curBar = colorbar;
        curBar.Ticks = linspace(0,255,6);
        curBar.TickLabels = arrayfun(@(x) sprintf('%0.4f',x),...
            min(totalFcn):0.2*(max(totalFcn)-min(totalFcn)):max(totalFcn),'un',0);
        cbh(i) = curBar;
    end
end

Link=linkprop(h, {'CameraUpVector','CameraPosition', 'CameraTarget',...
    'CameraViewAngle'});
end

