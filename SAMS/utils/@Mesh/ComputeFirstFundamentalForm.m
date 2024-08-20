function [E,F,G] = ComputeFirstFundamental(G,flatMesh)
% Old curvature methods based on Ronald Lui code. Not used but potentially
% interesting on their own
%Compute the Hopf differential of a mesh G with parametrization defined by
%G_uni. Experiemental method.
options.niter_averaging=4;
%options.type='combinatorial';
[x_Grad,~]=flatMesh.ComputeFunctionGradient(G.V(1,:)');
[y_Grad,~]=flatMesh.ComputeFunctionGradient(G.V(2,:)');
[z_Grad,~]=flatMesh.ComputeFunctionGradient(G.V(3,:)');
%smooth out so things aren't as ugly

% xCoord_u = G.PerformMeshSmoothing(xCoord_u,options);
% xCoord_v = G.PerformMeshSmoothing(xCoord_v,options);
% yCoord_u = G.PerformMeshSmoothing(yCoord_u,options);
% yCoord_v = G.PerformMeshSmoothing(yCoord_v,options);
% zCoord_u = G.PerformMeshSmoothing(zCoord_u,options);
% zCoord_v = G.PerformMeshSmoothing(zCoord_v,options);

xCoord_u = G.PerformMeshSmoothing(x_Grad(1,:),options);
xCoord_v = G.PerformMeshSmoothing(x_Grad(2,:),options);
yCoord_u = G.PerformMeshSmoothing(y_Grad(1,:),options);
yCoord_v = G.PerformMeshSmoothing(y_Grad(2,:),options);
zCoord_u = G.PerformMeshSmoothing(z_Grad(1,:)',options);
zCoord_v = G.PerformMeshSmoothing(z_Grad(2,:),options);

%differentiate again, pray
% [xCoord_uu,xCoord_uv]=trigradient(uv(:,1),uv(:,2),xCoord_u,G.F');
% [yCoord_uu,yCoord_uv]=trigradient(uv(:,1),uv(:,2),yCoord_u,G.F');
% [zCoord_uu,zCoord_uv]=trigradient(uv(:,1),uv(:,2),zCoord_u,G.F');
% [xCoord_vu,xCoord_vv]=trigradient(uv(:,1),uv(:,2),xCoord_v,G.F');
% [yCoord_vu,yCoord_vv]=trigradient(uv(:,1),uv(:,2),yCoord_v,G.F');
% [zCoord_vu,zCoord_vv]=trigradient(uv(:,1),uv(:,2),zCoord_v,G.F');


%time to compute the mu parameter. First, organize the parametrization
%derivatives

Ru = [xCoord_u';yCoord_u';zCoord_u'];
Rv = [xCoord_v';yCoord_v';zCoord_v'];

E = sum(Ru .* Ru);
F = sum(Ru .* Rv);
G = sum(Rv .* Rv);
end