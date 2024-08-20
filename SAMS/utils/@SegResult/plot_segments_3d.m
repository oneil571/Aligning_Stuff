function plot_segments_3d(SegResult, n, filePath)
% PLOT_SEGMENTS_3D - Plots and opt. saves figure of segments for n random meshes
% Adapted from original util/ViewBundleFunc.m

	disp('Plotting and saving representative collection of segments...');
	figure('Position', [10, 10, 800, 800]);
	set(gcf, 'ToolBar', 'none');
	h = zeros(size(n));
	rng('shuffle');
	hIdx = randi(length(SegResult.mesh), n, 1);

    if n/4 > 5
        nCol = ceil(n/4);
    else
        nCol = ceil(n/2);
    end
        nRow = ceil(n/nCol);
        

	for i=1:n
		m = SegResult.mesh{hIdx(i)};
	    color_data = m.segmentIdx;
	    h(i) = subplot(nRow, nCol, i);
	    m.draw(struct('FaceColor', 'interp', 'FaceVertexCData', color_data, 'CDataMapping','scaled', 'EdgeColor', 'none', 'FaceAlpha',1,'AmbientStrength',0.3,'SpecularStrength',0.0));
	    hold on;
	    colormap jet(256);
	    camlight('headlight');
	    camlight(180,0);
	    %title(m.Aux.Name);
	end

	link = linkprop(h, {'CameraUpVector', 'CameraPosition', 'CameraTarget', 'CameraViewAngle'});
	setappdata(gcf, 'StoreTheLink', link);

	savefig(gcf, filePath);

end