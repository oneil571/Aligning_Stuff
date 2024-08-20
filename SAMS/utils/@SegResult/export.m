function export(SegResult, cfg)
% EXPORT - Saves SegResult object and all data with write/plot methods

	outPath = cfg.out;

	if cfg.colorSegments
		SegResult.gen_segment_color();
	end

	%%% Save result object
	save(fullfile(outPath, 'result.mat'), 'SegResult');

	%%% Write mesh, segment, group PLYs
	SegResult.write_meshes(fullfile(outPath, 'mesh'));
	fclose('all');
	SegResult.write_segments(fullfile(outPath, 'segment'), ...
		1, cfg.dirCollate, cfg.colorSegments);
	fclose('all');
	SegResult.write_seg_all(fullfile(outPath, 'seg_all.off'), ...
		cfg.colorSegments);
	fclose('all');

	%%% Write CSV tables
	SegResult.write_mesh_table_csv(fullfile(outPath, 'mesh.csv'));
	SegResult.write_seg_table_csv(fullfile(outPath, 'segment.csv'));
	SegResult.write_seg_group_csvs(fullfile(outPath, 'seg_groups.csv'), fullfile(outPath, 'mesh_seg_groups.csv'));

	%%% Plot and save figures
	SegResult.plot_freq_dist(fullfile(outPath, 'freq_dist.eps'));
	SegResult.plot_segments_3d(cfg.numMeshDisplay, ...
		fullfile(outPath, 'segment.fig'));

	%%% Save result object
	save(fullfile(outPath, 'result.mat'), 'SegResult');

	disp('File export done!');

end