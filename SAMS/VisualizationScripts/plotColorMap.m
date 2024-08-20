close all
if ~exist([workingPath 'newMeshList.mat'])
    error('Reparametrized meshes do not exist, do not call until finished');
end
load([workingPath 'newMeshList.mat']);
load([workingPath 'FinalDists.mat']);
frechMean = find(sum(dists.^2) == min(sum(dists.^2)));
frechMean = frechMean(1);
frech=newMeshList{frechMean};
colorMap = frech.V;
colorMap(1,:) = (colorMap(1,:)-min(colorMap(1,:)))';
colorMap(2,:) = (colorMap(2,:)-min(colorMap(2,:)))';
colorMap(3,:) = (colorMap(3,:)-min(colorMap(3,:)))';
colorMap(1,:) = colorMap(1,:)/max(colorMap(1,:));
colorMap(2,:) = colorMap(2,:)/max(colorMap(2,:));
colorMap(3,:) = colorMap(3,:)/max(colorMap(3,:));

clear h
numRows = ceil(length(newMeshList)/visNumRow);
for i = 1:length(newMeshList)
h(i) = subplot(numRows,visNumRow,i);
newMeshList{i}.draw(struct('FaceColor', 'interp', 'FaceVertexCData', colorMap', 'CDataMapping', 'scaled', 'EdgeColor', 'none', 'FaceAlpha', 1, 'AmbientStrength',0.7,'SpecularStrength',0.0));
%scatter3(newMeshList{i}.V(1,:),newMeshList{i}.V(2,:),newMeshList{i}.V(3,:),40,colorMap','MarkerFaceColor','flat',...
    %'MarkerFaceAlpha',.5);
axis off
hold on
end
Link = linkprop(h, {'CameraUpVector', 'CameraPosition', 'CameraTarget', 'CameraViewAngle'});
setappdata(gcf, 'StoreTheLink', Link);
touch([workingPath 'MapFigures/']);
touch([workingPath 'MapFigures/FullMaps/']);
savefig([workingPath 'MapFigures/FullMaps/FullMaps.fig']);
saveas(gcf,[workingPath 'MapFigures/FullMaps/FullMaps.png']);