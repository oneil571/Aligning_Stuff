function write_segments(SegResult, dirPath, excludeEmpty, subdirByMesh, colorSegments)
% WRITE_SEGMENTS - Save segment meshes
	
	touch(dirPath);
	if subdirByMesh
		n = cellfun(@(x) x.Aux.name, SegResult.mesh, 'UniformOutput', 0);
		n = strrep(n, '.', '_');
		d = cellfun(@(x) fullfile(dirPath, x), n, 'UniformOutput', 0);
		cellfun(@touch, d);
	else
		d = repmat({dirPath}, 1, length(SegResult.mesh));
	end
		
	for i = 1:length(SegResult.mesh)
		disp(['Saving segments for mesh ' SegResult.mesh{i}.Aux.Name '...']);
		for j = 1:length(SegResult.mesh{i}.segment)
			if (excludeEmpty) && (size(SegResult.mesh{i}.segment{j}.F, 2) == 0) 
				continue
			end
			if colorSegments
				options.color = SegResult.mesh{i}.segment{j}.Aux.c;
			else
				options = struct;
			end
			SegResult.mesh{i}.segment{j}.Write(fullfile(d{i}, ...
				[SegResult.mesh{i}.Aux.Name '_seg' num2str(j) '.ply']), ...
				'ply', options);
			fclose('all');
		end
	end

end