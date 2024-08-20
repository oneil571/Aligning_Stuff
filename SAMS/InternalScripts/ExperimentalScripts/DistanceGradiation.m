param = 1%0:10:110;
dists = cell(length(param),1);

for i = 1:length(Names)
    disp(i)
	for j = 1:length(Names)
        normI = newMeshList{i}.ComputeNormal;
        normJ = newMeshList{j}.ComputeNormal;
        %normDis = mean(1-dot(normI,normJ));
        normDis = norm(normI-normJ)/newMeshList{i}.nV;
    	for k = 1:length(param)

            dists{k}(i,j) = 0*norm(newMeshList{i}.V - newMeshList{j}.V) + param(k) *normDis;
        end
    end
end
figure; clear h; hold on;

for i = 1:length(param)
    h(i) = subplot(ceil(length(param)/6),6,i);
    imagesc(dists{i})
end
Link = linkprop(h, {'CameraUpVector', 'CameraPosition', 'CameraTarget', 'CameraViewAngle'});