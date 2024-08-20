%% Set definitions

load([workingPath 'GPDists.mat']);
load([workingPath 'Names.mat']);
load([workingPath 'MappingData/MatchesPairs_Thresheld.mat']);
frechMean = find(min(sum(GPDists.^2))==sum(GPDists.^2));

samplesPath = [workingPath 'ProcessedMAT/'];
meshes = cell(length(Names),1);

disp('Loading meshes...')
for i = 1:length(Names)
    load([samplesPath Names{i} '.mat']);
    meshes{i} = G;
end
meshList = meshes;
options.FeatureType = 'ConfMax';
options.NumDensityPnts = 100;
options.AngleIncrement = 0.01;
options.NumFeatureMatch = 4;
options.GaussMinMatch = 'off';

TextureCoords1 = cell(length(Names),1);
TextureCoords2 = cell(length(Names),1);
TextureCoordsRev1 = cell(length(Names),1);
TextureCoordsRev2 = cell(length(Names),1);
maps = cell(length(Names),1);
translation = cell(length(Names),1);


%% Compute initial alignments based on landmarks
progressbar
for i = 1:length(Names)
    progressbar(i/length(Names));
    if i == frechMean
        continue;
    end
    if isempty(matchesPairs{i})
        disp('No landmarks matched, continuing');
        continue;
    end
    options.GPMatches = matchesPairs{i};
    curMatchedLmks = matchesPairs{i};
    %Gather points
    ptCloud_1 = meshes{i}.V(:,curMatchedLmks(:,1));
    ptCloud_2 = meshes{frechMean}.V(:,curMatchedLmks(:,2));

    %centralize and normalize
    ptCloud_1 = ptCloud_1 - repmat(mean(ptCloud_1')',1,size(ptCloud_1,2));
    ptCloud_2 = ptCloud_2 - repmat(mean(ptCloud_2')',1,size(ptCloud_2,2));
    ptCloud_1 = ptCloud_1/norm(ptCloud_1,'fro');
    ptCloud_2 = ptCloud_2/norm(ptCloud_2,'fro');

    [U,~,V] = svd(ptCloud_1*(ptCloud_2'));
    for q = 1:meshes{i}.nV
        meshes{i}.V(:,q) = V*U'*meshes{i}.V(:,q);
    end
end
close all
R = cell(length(Names),1);
%% Compute registrations
disp('Computing registrations');
isoParam = cell(513,1);
i=frechMean;

[meshList{i}.Aux.V0, meshList{i}.Aux.inds_bF] = computeTutte(meshList{i}.F',meshList{i}.V',false); % global scale to be as ismetric as possible
% constrain the centroid
meshList{i}.Aux.eq_lhs_centroid = kron(eye(2),sparse(1,meshList{i}.Aux.inds_bF,1,1,size(meshList{i}.V,2)));
meshList{i}.Aux.eq_rhs_centroid = [0;0];
% setup optimization problem
optimProblemParam = OptimProblemIsoDist(meshList{i}.V', meshList{i}.F', meshList{i}.Aux.eq_lhs_centroid, meshList{i}.Aux.eq_rhs_centroid, meshList{i}.Aux.V0);
% setup solver
solverParam = OptimSolverAcclQuadProx('QP', optimProblemParam, false, true, true);
% solve
logParam = solverParam.solveTol(TolX, TolFun ,num_iter);
% return isometric parameterization
frechUni = Mesh('VF',[solverParam.X';zeros(1,meshList{i}.nV)],meshList{i}.F);
    

progressbar
for i = 403:length(Names)
    %progressbar(i/length(Names));
    if i == frechMean
        TextureCoords1{i} = frechUni.V(1:2,:);
        TextureCoords2{i} = TextureCoords1{i};
        TextureCoordsRev1{i} = TextureCoords1{i};
        TextureCoordsRev2{i} = TextureCoords1{i};
        continue;
    end
    [meshList{i}.Aux.V0, meshList{i}.Aux.inds_bF] = computeTutte(meshList{i}.F',meshList{i}.V',false); % global scale to be as ismetric as possible
    % constrain the centroid
    meshList{i}.Aux.eq_lhs_centroid = kron(eye(2),sparse(1,meshList{i}.Aux.inds_bF,1,1,size(meshList{i}.V,2)));
    meshList{i}.Aux.eq_rhs_centroid = [0;0];
    % setup optimization problem
    optimProblemParam = OptimProblemIsoDist(meshList{i}.V', meshList{i}.F', meshList{i}.Aux.eq_lhs_centroid, meshList{i}.Aux.eq_rhs_centroid, meshList{i}.Aux.V0);
    % setup solver
    solverParam = OptimSolverAcclQuadProx('QP', optimProblemParam, false, true, true);
    % solve
    logParam = solverParam.solveTol(TolX, TolFun ,num_iter);
    % return isometric parameterization
    curUni = Mesh('VF',[solverParam.X';zeros(1,meshList{i}.nV)],meshList{i}.F);
    curPairs = matchesPairs{i};
   

    % global alignment -- apply global rigid tranformation to {1} based on correspondences
    [T.R, T.t] = findBestRigidMotion(curUni.V(1:2,curPairs(:,1))',frechUni.V(1:2,curPairs(:,2))',1);
    % compute linear system (corresponding to landmark matches)
    V = (curUni.V(1:2,:)'+T.t)*T.R';
    if det(T.R)>0
        F =  curUni.F';
    else
        F =  curUni.F([1 3 2],:)';
    end
    [eq_lhs,eq_rhs] = indCoordsToLinearSystem(curUni.V(1:2,:)',...
        curPairs(:,1), ...
        frechUni.V(1:2,curPairs(:,2))');
    % setup optimization problem
    num_iter = 500;
    TolX = 1e-10;
    TolFun = 1e-6;
    optimProblem = OptimProblemIsoDist(V, F, eq_lhs, eq_rhs, [], 5);
    % setup solver
    solverInterp = OptimSolverAcclQuadProx('QP', optimProblem, false, true, true);
    % solve
    logInterp = solverInterp.solveTol(TolX, TolFun ,num_iter);
    % store registered parameterizations
    
    
    %rslt_GP = meshes{i}.InterpolateCorrespondences(meshes{frechMean},options);
    %rslt_cP = meshes{i}.ComputeContinuousProcrustes(meshes{frechMean},options);
    
    %if rslt_GP.cPdist < rslt_cP.cPdist
        TextureCoords1{i} = solverInterp.X';
        TextureCoords2{i} = frechUni.V;
        %maps{i} = rslt_GP.cPmap;
        %translation{i} = rslt_GP.translation;
        %R{i} = rslt_GP.orthogonal;


end
close all

save([workingPath 'TextureCoords1.mat','TextureCoords1']);
save([workingPath 'TextureCoords2.mat','TextureCoords2']);
% figure;
% for i = 1:length(Names)
% h(i) = subplot(5,11,i);
% meshes{i}.draw;
% hold on
% end

