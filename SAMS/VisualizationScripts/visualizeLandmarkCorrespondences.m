try
    load([workingPath 'MappingData/MatchesPairs_Thresheld.mat']);
catch
    error('Please compute landmark correspondences before running this method');
end
disp('Visualizing all distilled correspondences');
disp('Press any key to switch chunks');
load([workingPath 'Names.mat']);
visMeshList = cell(length(Names),1);
for i = 1:length(Names)
    load([workingPath 'ProcessedMAT/' Names{i} '.mat']);
    visMeshList{i} = G;
end
numChunk = ceil((length(Names)-1)/visNumRow);
load([workingPath 'GPDists.mat']);
frechMean = find(sum(GPDists.^2)==min(sum(GPDists.^2)));
frechMesh = visMeshList{frechMean};
visMeshList(frechMean) = [];
visMatchesPairs = matchesPairs([1:frechMean-1 frechMean+1:length(Names)]);
for i = 1:numChunk
    startInd =1+(i-1)*visNumRow;
    figure
    clear h
    for j = startInd:min(startInd+visNumRow-1,length(Names)-1)
        h(1,j-startInd+1) = subplot(2,visNumRow,j-startInd+1);
        visMeshList{j}.draw; hold on;
        curV = visMeshList{j}.V;
        curCube = distinguishable_colors(size(visMatchesPairs{j},1));
        scatter3(curV(1,visMatchesPairs{j}(:,1)),curV(2,visMatchesPairs{j}(:,1)),...
            curV(3,visMatchesPairs{j}(:,1)),100,curCube,'filled');
        
        h(2,j-startInd+1) = subplot(2,visNumRow,j-startInd+1+visNumRow);
        frechMesh.draw; hold on;
        curV = frechMesh.V;
        scatter3(curV(1,visMatchesPairs{j}(:,2)),curV(2,visMatchesPairs{j}(:,2)),...
            curV(3,visMatchesPairs{j}(:,2)),100,curCube,'filled');
    end
    h = reshape(h,length(h(:)),1);
    Link = linkprop(h, {'CameraUpVector', 'CameraPosition', 'CameraTarget', 'CameraViewAngle'});
    setappdata(gcf, 'StoreTheLink', Link);
    pause()
    close all
end
