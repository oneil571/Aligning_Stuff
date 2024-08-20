function write_seg_detail_csv(SegResult, filePath)
% WRITE_SEG_DETAIL_CSV - Write CSV with segment table

	disp('Writing CSV table of segment details...');
	d = SegResult.data;
	meshName    = {};
	segNum      = [];
	segNFace    = [];
	segNVert    = [];
	segArea     = [];
	segPercArea = [];
    segDNE      = [];
    segPercDNE  = [];

	for i = 1:length(SegResult.mesh)
		for j = 1:length(SegResult.mesh{i}.segment)
			meshName = [meshName {d.meshName{i}}];
			segNum = [segNum j];
			segNFace = [segNFace size(SegResult.mesh{i}.segment{j}.F, 2)];
			segNVert = [segNVert size(SegResult.mesh{i}.segment{j}.V, 2)];
			a = d.segmentArea{i}(j);
            b = d.segmentDNE{i}(j);
			segArea = [segArea a];
			segPercArea = [segPercArea a/d.segmentAreaTotal(i)]; 
            segDNE = [segDNE b];
            segPercDNE = [segPercDNE b/d.segmentDNETotal(i)];
		end
	end

	varNames = {'mesh_name', 'segment_id', 'n_face', 'n_vert', 'area', ...
		'proportional_area','dne','proportional_dne'};
	dataTable = table(meshName', segNum', segNFace', segNVert', segArea', ...
		segPercArea', segDNE',segPercDNE', 'VariableNames', varNames);
	writetable(dataTable, filePath);

end