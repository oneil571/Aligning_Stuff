function write_meshes(SegResult, dirPath)
% WRITE_MESHES - Save mesh files
	
	touch(dirPath);

	for i = 1:length(SegResult.mesh)
		disp(['Saving mesh ' SegResult.mesh{i}.Aux.Name ' as OFF file...']);
		SegResult.mesh{i}.Write(fullfile(dirPath, ...
			[SegResult.mesh{i}.Aux.Name '.ply']), 'ply', struct);
		fclose('all');
	end
	
end