function [newG,Aux] = ComputeAuxiliaryInformation(G,options)
%Wrapper function to compute relevant quantitative auxilliary quantities. Needed in
%order to do in parallel fashion. This will set all auxiliary data except
%for the name.

num_iter = 500;
TolX = 1e-10;
TolFun = 1e-6;
%% Get aux field if already exists, should only be if observer landmarks were computed
if isfield(G,'Aux')
    Aux = G.Aux;
else
    Aux = {};
end

%centralize
[Aux.Area,Aux.Center] = G.Centralize('ScaleArea');

%uniformizae and compute conformal factor
switch options.isDisc
    case 1
        [Aux.UniformizationV,~,Aux.VertArea] = G.ComputeMidEdgeUniformization();
        
        [G.Aux.V0, G.Aux.inds_bF] = computeTutte(G.F',G.V',false); % global scale to be as ismetric as possible
% constrain the centroid
G.Aux.eq_lhs_centroid = kron(eye(2),sparse(1,G.Aux.inds_bF,1,1,size(G.V,2)));
G.Aux.eq_rhs_centroid = [0;0];
% setup optimization problem
optimProblemParam = OptimProblemIsoDist(G.V', G.F', G.Aux.eq_lhs_centroid, G.Aux.eq_rhs_centroid, G.Aux.V0);
% setup solver
solverParam = OptimSolverAcclQuadProx('QP', optimProblemParam, false, true, true);
logParam = solverParam.solveTol(TolX, TolFun ,num_iter);
% return isometric parameterization
Aux.IsoUniV = [solverParam.X';zeros(1,G.nV)];
    case 0
        [~,TriAreas] = G.ComputeSurfaceArea;
        Aux.VertArea = (TriAreas'*G.F2V)/3;
        Aux.UniformizationV = spherical_conformal_map(G.V',G.F')';
end

disp('Computing curvatures and Laplacian');
Aux = G.ComputeCurvatureFeatures(Aux,options);
Aux.LB = G.ComputeCotanLaplacian;
G.Aux = Aux;

disp('Computing curvature maxima');
minds = [Aux.GaussMaxInds;Aux.MeanMinInds;Aux.ConfMaxInds;Aux.DNEMaxInds];
minds = unique(minds);

%numDensityPts = getoptions(options,'NumDensityPts',100);
disp('Computing GP Landmarks')
[G.Aux.GPLmkInds,~] = G.GetGPLmk(options.numGPLmks);
%disp('Computing FPS Landmarks')
%G.Aux.DensityPnts = G.GeodesicFarthestPointSampling(numDensityPts,minds);
newG = G;
Aux = G.Aux;
disp('All mesh information computed!');
end