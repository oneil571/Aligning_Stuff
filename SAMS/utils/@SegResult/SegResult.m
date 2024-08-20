classdef SegResult < handle
% SEGRESULT - Class storing results from surface region spectral clustering
% This class stores surface meshes, and generates and stores mesh surface region 
% segments given a concatenated column vector of segment ID values (1:n) per 
% vertex for all surface meshes. Methods calculate secondary data and export 
% surface meshes and secondary data. 
%
% Constructor syntax: obj = SegResult(flatSamples, kIdx, vIdxCumSum, cfg)
% 
% Constructor inputs: 
%	flatSamples -  
%	kIdx        -
%	vIdxCumSum  -
%	cfg         -
% 
% Attributes:
%	cfg  -
%	data -
%	mesh -
%	R    -



	properties
	% Documentation needed
		data = struct;
		mesh = {};
		R = {};
	end

	methods

		function obj = SegResult(flatSamples, kIdx, vIdxCumSum)
		% Class constructor
		
			vIdxCumSum = [0; vIdxCumSum];
			for i = 1:length(flatSamples)
				segIdx = kIdx(vIdxCumSum(i)+1:vIdxCumSum(i+1));
			    tmp = flatSamples{i}; 
			    obj.mesh{i} = SegMesh(tmp);
			    obj.mesh{i}.segmentIdx = segIdx;
			    obj.mesh{i}.segment = cell(max(segIdx), 1);
			    V = obj.mesh{i}.V';
			    F = obj.mesh{i}.F';
			    for j = 1:size(obj.mesh{i}.segment, 1)
			    	obj.mesh{i}.segment{j} = Mesh();
			        obj.mesh{i}.segment{j}.V = V(segIdx == j, :)';
			        vIdx = 1:length(V);
			        segVertIdxOrig = vIdx(segIdx == j);
			        segVertIdxNew = 1:length(segVertIdxOrig);
			        segFaceIdxOrig = F(all(ismember(F, segVertIdxOrig), 2),:);
			        segFaceIdxNew = changem(segFaceIdxOrig, segVertIdxNew, ...
			        	segVertIdxOrig);
			        obj.mesh{i}.segment{j}.F = segFaceIdxNew';
			    end
            end
			
		end

	end
end