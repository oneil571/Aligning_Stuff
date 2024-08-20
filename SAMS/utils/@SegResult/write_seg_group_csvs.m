function write_seg_group_csvs(SegResult, segGroupPath, meshSegGroupPath)
% WRITE_SEG_GROUP_CSVS - Saves SegResult segGroup object tables

	disp('Writing CSV tables of segment distribution groups...');
	g = SegResult.data.segGroups;

	% Segment group csv
	names = keys(g)';
	lens = cellfun(@(x) length(strsplit(x))-1, names, 'UniformOutput', false);
	mesh_num = values(g)';
	segGroupVarNames = {'segment_group', 'segments_number', 'mesh_number'};
	segGroupTable = table(names, lens, mesh_num, 'VariableNames', segGroupVarNames);
	writetable(segGroupTable, segGroupPath);

	meshSegGroups = cell(length(SegResult.mesh), 1);
	meshSegGroupIdxs = cell(length(SegResult.mesh), 1);
	% Meshes by segment group csv
	for i = 1:length(SegResult.mesh)
		meshSegGroups{i} = SegResult.mesh{i}.segGroup;
		meshSegGroupIdxs{i} = find(strcmp(names,meshSegGroups{i}));
	end
	meshSegGroupVarNames = {'mesh_name', 'segment_group', 'segment_group_id'};
	meshSegGroupTable = table(SegResult.data.meshName, meshSegGroups, meshSegGroupIdxs, 'VariableNames', meshSegGroupVarNames);
	writetable(meshSegGroupTable, meshSegGroupPath);

end