function h = DrawTextureMap(uv,xyz,I,scale)
%% Written by Shahar Kovalsky

if ~exist('scale','var')
    scale = 1;
end
uv_range = [0 1 0 1]*scale;
[u_grid, v_grid] = meshgrid(linspace(uv_range(1),uv_range(2),size(I,1)), linspace(uv_range(3),uv_range(4),size(I,2)));

F = cell(3,1);
xyz_intrp = cell(3,1);
for jj = 1:3
    F{jj} = scatteredInterpolant(uv,xyz(:,jj),'linear','none');
    xyz_intrp{jj} = F{jj}(u_grid, v_grid);
end

h = warp(xyz_intrp{1}, xyz_intrp{2}, xyz_intrp{3}, I);
axis off
axis tight
axis equal
         set(gca, 'CameraViewAngle', 10.0659);
        camlight('headlight');
        camlight(180,0);

