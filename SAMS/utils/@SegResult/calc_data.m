function calc_data(SegResult)
% Calculates mesh/segment properties and populates SegResult.data fields

	r = SegResult;

	data = struct;
	data.meshN            = length(SegResult.mesh);
	data.meshArea         = zeros(data.meshN, 1);
    data.meshDNE         = zeros(data.meshN, 1);
	data.meshBorderArea   = zeros(data.meshN, 1);
	data.segmentTotalN    = zeros(data.meshN, 1);
	data.segmentNonZeroN  = zeros(data.meshN, 1);
	data.segmentAreaTotal = zeros(data.meshN, 1);
	data.segmentArea      = cell(data.meshN, 1);
    data.segmentDNETotal = zeros(data.meshN, 1);
	data.segmentDNE     = cell(data.meshN, 1);
	data.meshName		  = cellfun(@(x) x.Aux.Name, ...
		SegResult.mesh, 'UniformOutput', 0)';

	for i = 1:data.meshN		
		data.meshArea(i) = r.mesh{i}.Aux.Area;
        rslt = r.mesh{i}.ariaDNE; data.meshDNE(i) = rslt.dne;
		data.segmentTotalN(i) = length(r.mesh{i}.segment);
		data.segmentNonZeroN(i) = ...
			nnz(cellfun(@(x) size(x.F, 2), r.mesh{1}.segment));
		segAreaArray = zeros(data.segmentTotalN(i), 1);
        segDNEArray = zeros(data.segmentTotalN(i),1);
        curIdxs = r.mesh{i}.segmentIdx;
		for j = 1:data.segmentTotalN(i)
			segAreaArray(j) = data.meshArea(i)*surface_area(r.mesh{i}.segment{j}.V, ...
				r.mesh{i}.segment{j}.F);
            curSegInds = find(curIdxs == j);
            segDNEArray(j) = sum(rslt.localDNE(curSegInds));
		end
		data.segmentArea{i} = segAreaArray;
        data.segmentDNE{i} = segDNEArray;
		data.segmentAreaTotal(i) = sum(segAreaArray); 
        data.segmentDNETotal(i) = sum(segDNEArray);
	end
	data.meshBorderArea = data.meshArea - data.segmentAreaTotal;
	SegResult.data = data;
	
	SegResult.collate_seg_groups();

	function A = surface_area(V, F)
		% Expects dim x nV and dim x nF matrices
		pArea = zeros(size(F, 2), 1);
		for k = 1:length(pArea)
			p = V(:, F(:, k));
			pArea(k) = triangle_area(p(:, 1), p(:, 2), p(:, 3));
		end
		A = sum(pArea);
	end

	function A = triangle_area(a, b, c)
		A = 0.5 * sqrt(sum(abs(cross(b-a, c-a)).^2));
	end

end