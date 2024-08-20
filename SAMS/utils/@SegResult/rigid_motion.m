function rigid_motion(SegResult)
% RIGID_MOTION - flatSamples x flatSamples cell array of rotation matrices

	load(fullfile(SegResult.cfg.path.cpdMST, 'r_mst.mat'));
	SegResult.R = r;

end