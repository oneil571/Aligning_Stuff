close all
addpath(genpath('./Mapping/external/AQP_toolbox/AcceleratedQuadraticProxy/'));

load([workingPath 'Flags.mat']);
if ~exist([workingPath 'newMeshList.mat'])
    error('Reparametrized meshes do not exist, do not call until finished');
elseif ~Flags('isDisc')
    error('This method only works for surfaces with disk topology, aborting...')
end
    
load([workingPath 'newMeshList.mat']);
for i = 1:length(newMeshList)
    newMeshList{i}.V = newMeshList{i}.V-mean(newMeshList{i}.V,2);
    newMeshList{i}.V = newMeshList{i}.V/norm(newMeshList{i}.V,'fro');
end
load([workingPath 'FinalDists.mat']);
frechMean = find(sum(dists.^2) == min(sum(dists.^2)));
frechMean = frechMean(1);
frechMesh = newMeshList{frechMean};

if ~isfield(frechMesh.Aux,'UniformizationV')
    disp('Template mesh does not have texture coordinates')
    disp('Computing texture coordinates and other features')
    options.isDisc = 1; options.numGPLmks=200;
    frechMesh.ComputeAuxiliaryInformation(options);
end
frechMesh.V = frechMesh.V-mean(frechMesh.V,2);
frechMesh.V = frechMesh.V/norm(frechMesh.V,'fro');
TextureCoords = frechMesh.Aux.IsoUniV(1:2,:)';
TextureCoords = (TextureCoords-min(min(TextureCoords))+0.01);
TextureCoords = TextureCoords/(max(max(TextureCoords))+0.01);
%%
clear h
%numRows = ceil(length(newMeshList)/visNumRow);
uv_grid=imread('uv_grid.jpg');
uv_grid = uv_grid(:,end:-1:1,:);
negFlag = 1;
disp('Drawing texture maps...')
figure; hold on; axis off
%progressbar
%F = newMeshList{1}.F;
%if negFlag == 1
%    F = [F(1,:);F(3,:);F(2,:)]';
%end
for i = 1:10%112:116%length(newMeshList)
    %progressbar(i/length(newMeshList));
    %disp([num2str(i) '/' num2str(length(newMeshList))])
    h(i) = subplot(2,5,i);%numRows,visNumRow,i);
    DrawTextureMap(TextureCoords,newMeshList{i}.V',uv_grid);
    %patchtexture(meshList{i}.F',meshList{i}.V',meshList{i}.F',TextureCoords(mapList{i,frechMean},:),uv_grid);
end
Link = linkprop(h, {'CameraUpVector', 'CameraPosition', 'CameraTarget', 'CameraViewAngle'});
setappdata(gcf, 'StoreTheLink', Link);
touch([workingPath 'MapFigures/']);
touch([workingPath 'MapFigures/FullMaps/']);
%disp('Saving Texture Maps...')
%savefig(h,[workingPath 'MapFigures/FullMaps/TextureMaps.fig']);
%saveas(gcf,[workingPath 'MapFigures/FullMaps/TextureMaps.png']);
disp('Done!')