function collate_seg_groups(SegResult)
% Collates segmentation distribution groups for result and per mesh

	r = SegResult;
	segGroups = containers.Map;
	for i = 1:length(r.mesh)
		% Determine mesh segment distribution group (i.e., "1 3 5")
		meshSegGroup = '';
		for j = 1:length(r.mesh{i}.segment)
			if length(r.mesh{i}.segment{j}.V) > 0
				meshSegGroup = [meshSegGroup, num2str(j), ' '];
			end
		end
		r.mesh{i}.segGroup = meshSegGroup;

		% Determine whether group is known or must be created
		% and assort mesh to this group
		if ismember(meshSegGroup, keys(segGroups))
			segGroups(meshSegGroup) = segGroups(meshSegGroup) + 1;
		else
			segGroups(meshSegGroup) = 1;
		end
	end
	r.data.segGroups = segGroups;

end