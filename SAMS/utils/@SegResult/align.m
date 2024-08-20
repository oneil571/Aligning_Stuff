function align(SegResult, idx)
% ALIGN - Aligns all SegResult meshes and segments to idx mesh

	SegResult.rigid_motion();
	for i = 1:length(SegResult.mesh)
		SegResult.mesh{i}.V = SegResult.R{i, idx} * SegResult.mesh{i}.V;
		if det(SegResult.R{i, idx}) < 0
			SegResult.mesh{i}.F = flipud(SegResult.mesh{i}.F);
		else
			SegResult.mesh{i}.F = SegResult.mesh{i}.F;
		end
		for j = 1:length(SegResult.mesh{i}.segment)
			SegResult.mesh{i}.segment{j}.V = SegResult.R{i, idx} * ...
				SegResult.mesh{i}.segment{j}.V;
			if det(SegResult.R{i, idx}) < 0
				SegResult.mesh{i}.segment{j}.F = ...
					flipud(SegResult.mesh{i}.segment{j}.F);
			else
				SegResult.mesh{i}.segment{j}.F = SegResult.mesh{i}.segment{j}.F;
			end
		end
	end

end