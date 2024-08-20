function write_seg_all(SegResult, filePath, colorSegments)
% Creates and saves a mesh combining all segments from all meshes for comparison

	disp('Writing whole segment sample OFF...');
	newMesh = SegResult.mesh;
	
	nMesh = length(newMesh);
	nRow = floor(sqrt(nMesh));
	nCol = ceil(nMesh/nRow);
	xMax = max(cellfun(@(m) max(m.V(1, :)) - min(m.V(1, :)), newMesh));
	yMax = max(cellfun(@(m) max(m.V(2, :)) - min(m.V(2, :)), newMesh));
	xStep = xMax * 1.25;
	yStep = yMax * 1.25;

	% Locations of meshes
	xLocs = repmat(0:xStep:xStep*(nCol-1), 1, nRow);
	yLocs = reshape(repmat(0:yStep:yStep*(nRow-1), nCol, 1), [1, nRow*nCol]);
	zLocs = zeros(1, nRow * nCol);
	locs = vertcat(xLocs, yLocs, zLocs);

	groupV = [];
	groupF = [];
	groupC = [];
	for i = 1:length(newMesh)
		vectCenter = mean(newMesh{i}.V, 2);
		vectLoc = locs(:, i) - vectCenter;
		for j = 1:length(newMesh{i}.segment)
			segLoc = repmat(vectLoc, 1, size(newMesh{i}.segment{j}.V, 2));
			movedV = newMesh{i}.segment{j}.V + segLoc;
			groupF = [groupF newMesh{i}.segment{j}.F+length(groupV)];
			groupV = [groupV movedV];
			groupC = [groupC newMesh{i}.segment{j}.Aux.c];
		end
	end

	meshGroup = Mesh('VF', groupV, groupF);
	if colorSegments
		options.color = groupC;
	else
		options = struct;
    end
    options.pointCloud=0;
	meshGroup.Write(fullfile(filePath), 'off', options);

end

