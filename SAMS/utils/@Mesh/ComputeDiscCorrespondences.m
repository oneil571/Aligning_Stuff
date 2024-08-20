function [rslt] = ComputeContinuousProcrustesStable_FixedRotation(GM,GN,options)
%COMPUTECONTINUOUSPROCRUSTES: Compute cP distance between GM and GN. Output
%contians
%   rslt.Gname1:            name of the first mesh
%   rslt.Gname2:            name of the second mesh
%   rslt.cPdist:            continuous Procrustes distance
%   rslt.cPmap:             optimal map generating cP distance
%   rslt.invcPmap:          inverse of rslt.cPmap
%   rslt.TextureCoords1:    texture coordinates for the first mesh
%                           (deformed)
%   rslt.TextureCoords2:    textrue coordinates for the second mesh
%                           (not deformed)
%   rslt.ref:               =0 if cP map is orientation-preserving
%                           =1 if cP map is orientation-reversing
%
%   Tingran Gao, trgao10@math.duke.edu
%   Modified by Robert Ravier to include ability to match pre determined
%   user landmarks
%


if nargin<3
    options.FeatureType = 'ConfMax';
    options.NumGPLmks = 250;    %Should be big enough to make sure all features captured plus some cover of interior
    options.MaxDistTol = 8;    %This will be replaced with weighted adjacency in future
    options.NumDensityPnts = 100;
end
ProgressBar = getoptions(options,'ProgressBar','on');

%%% feature type for matching
FeatureType = getoptions(options,'FeatureType','ConfMax');
NumGPLmks = getoptions(options,'NumDensityPnts',250);
MaxDistTol = getoptions(options,'MaxDistTol',7);
NumDensityPnts = getoptions(options,'NumDensityPts',100);
switch FeatureType
    case 'ADMax'
        FeaturesM = GM.Aux.ADMaxInds;
        FeaturesN = GN.Aux.ADMaxInds;
    case 'GaussMax'
        FeaturesM = GM.Aux.GaussMaxInds;
        FeaturesN = GN.Aux.GaussMaxInds;
    case 'ConfMax'
        FeaturesM = GM.Aux.ConfMaxInds;
        FeaturesN = GN.Aux.ConfMaxInds;
    case 'Landmarks'
        FeaturesM = GM.Aux.Landmarks;
        FeaturesN = GN.Aux.Landmarks;
end

if length(FeaturesM) < 3 || length(FeaturesN) < 3
    disp('WARNING: Not enough features for matching. Computing candidate GP features');
    FeaturesM = GM.GetGPLmk(10);
    FeaturesN = GN.GetGPLmk(10);
end
FeaturesMCoords = compl(GM.Aux.UniformizationV(:,FeaturesM));
FeaturesNCoords = compl(GN.Aux.UniformizationV(:,FeaturesN));
map_12 = knnsearch(GN.V(:,FeaturesN)',GM.V(:,FeaturesM)');
map_21 = knnsearch(GM.V(:,FeaturesM)',GN.V(:,FeaturesN)');
if ~isfield(GM.Aux,'GPLmkInds')
    GM.Aux.GPLmkInds = GM.GetGPLmk(NumGPLmks);
end
if ~isfield(GN.Aux,'GPLmkInds')
    GN.Aux.GPLmkInds = GN.GetGPLmk(NumGPLmks);
end
%Coarse mapping: Nearest Neighbors
proc12 = knnsearch(GN.V(:,GN.Aux.GPLmkInds)',GM.V(:,GM.Aux.GPLmkInds)');
proc21 = knnsearch(GM.V(:,GM.Aux.GPLmkInds)',GN.V(:,GN.Aux.GPLmkInds)');
GP_M = GM.Aux.GPLmkInds;
GP_N = GN.Aux.GPLmkInds;
featureProj1 = knnsearch(GM.V(:,GM.Aux.GPLmkInds)',GM.V(:,FeaturesM)');
featureProj2 = knnsearch(GN.V(:,GN.Aux.GPLmkInds)',GN.V(:,FeaturesN)');
lmk12 = proc12(featureProj1);
lmk21 = proc21(featureProj2);
while 1 > 0
    TPS_FEATURESM_INDS = []; TPS_FEATURESN_INDS = [];
    for j = 1:length(map_12)
        if map_21(map_12(j)) == j
            curLmk12 = lmk12(j);
            curLmk21 = lmk21(map_12(j));
            [dist12,~,~] = graphshortestpath(GN.A,GP_N(curLmk12),GP_N(featureProj2(map_12(j))));
            [dist21,~,~] = graphshortestpath(GM.A,GP_M(curLmk21),GP_M(featureProj1(j)));
            if dist12 < MaxDistTol && dist21 < MaxDistTol
                TPS_FEATURESM_INDS = [TPS_FEATURESM_INDS FeaturesM(j)];
                TPS_FEATURESN_INDS = [TPS_FEATURESN_INDS FeaturesN(map_12(j))];
            end
        end
    end
    if length(TPS_FEATURESM_INDS) < 3
        disp(['Error: Not enough allowable correspondences. Expanding radius to '...
            num2str(MaxDistTol+1);]);
        MaxDistTol = MaxDistTol+1;
    else
        break;
    end
end
%%% check for NaN's in the uniformization of GM
sourceInds = GM.Aux.GPLmkInds(1:NumDensityPnts);
source = compl(GM.Aux.UniformizationV(:,sourceInds));
delInds = isnan(source);
source(delInds) = [];
sourceInds(delInds) = [];
VorArea = GM.ComputeVoronoiArea(sourceInds);
%%% check for NaN's in the uniformization of GN
targetInds = GN.Aux.GPLmkInds(1:NumDensityPnts);
target = compl(GN.Aux.UniformizationV(:,targetInds));
delInds = isnan(target);
target(delInds) = [];
targetInds(delInds) = [];


TextureCoords2 = GN.Aux.UniformizationV(1:2,:);
TextureCoords2(:,isnan(compl(TextureCoords2))) = ones(2,sum(isnan(compl(TextureCoords2))));
TextureCoords1 = GM.Aux.UniformizationV(1:2,:);

TPS_FEATURESM = DISCtoPLANE(TextureCoords1(:,TPS_FEATURESM_INDS)','d2p');
TPS_FEATURESN = DISCtoPLANE(TextureCoords2(:,TPS_FEATURESN_INDS)','d2p');

[TPS_FEATURESM,ia] = unique(TPS_FEATURESM,'rows');
TPS_FEATURESN = TPS_FEATURESN(ia,:);
[TPS_FEATURESN,ia] = unique(TPS_FEATURESN,'rows');
TPS_FEATURESM = TPS_FEATURESM(ia,:);
orig = TextureCoords1;
if (length(TPS_FEATURESM)>3) % TPS (Thin Plate Spline)
    tP = DISCtoPLANE(TextureCoords1','d2p');
    [ftps] = TEETH_calc_tps(TPS_FEATURESM,TPS_FEATURESN-TPS_FEATURESM);
    pt = tP + TEETH_eval_tps(ftps,tP);
    TextureCoords1 = DISCtoPLANE(pt,'p2d')';
elseif (length(TPS_FEATURESM)==3) % Affine Transformation
    tP = DISCtoPLANE(TextureCoords1','d2p');
    [A,b] = PlanarThreePtsDeform(TPS_FEATURESM,TPS_FEATURESN);
    pt = [A,b]*[tP';ones(1,size(tP,1))];
    TextureCoords1 = DISCtoPLANE(pt','p2d')';
end

%%% linearly interpolate texture coordinates for boundary points
if ~isfield(GM,'BV')
    [GM.BV,GM.BE] = GM.FindBoundaries();
end
size(orig)
THETA = cart2pol(orig(1,GM.BV),orig(2,GM.BV));
regIdx = find(~isnan(compl(TextureCoords1(:,GM.BV))));
nanIdx = find(isnan(compl(TextureCoords1(:,GM.BV))));
newTHETA = cart2pol(TextureCoords1(1,GM.BV(regIdx)),TextureCoords1(2,GM.BV(regIdx)));
%%%%%%%% THETA(regIdx) ---> newTHETA
%%%%%%%% THETA(nanIdx) ---> ??
interpBVTextureCoords = interp1(THETA(regIdx),newTHETA,THETA(nanIdx),'spline');
[X,Y] = pol2cart(interpBVTextureCoords,1);
TextureCoords1(:,GM.BV(nanIdx)) = [X;Y];
TextureCoords1(:,isnan(compl(TextureCoords1))) = orig(:,isnan(compl(TextureCoords1)));
%TextureCoords2_kdtree = kdtree_build(TextureCoords2');
cPmap = knnsearch(TextureCoords2', TextureCoords1');
%%% construct inverse map
%TextureCoords1_kdtree = kdtree_build(TextureCoords1');
invcPmap = knnsearch(TextureCoords1', TextureCoords2');

%%%%%
nBV = setdiff(1:GM.nV,GM.BV);

%%%%% linear interpolate images under the map when evaluating the cP functional
TR = triangulation(GN.F',TextureCoords2');
ti = pointLocation(TR,TextureCoords1(:,nBV)');
nBV(isnan(ti)) = [];
ti(isnan(ti)) = [];
BC = cartesianToBarycentric(TR,ti,TextureCoords1(:,nBV)');

imagPts = zeros(3,length(nBV));
for j=1:length(nBV)
    imagPts(:,j) = GN.V(:,GN.F(:,ti(j)))*BC(j,:)';
end

[cPdist,R,T] = MapToDist(GM.V(:,nBV),imagPts,1:length(nBV),GM.Aux.VertArea(nBV));
% cPdist = sqrt(sum((GM.V(:,nBV)-imagPts).^2)*GM.Aux.VertArea(nBV)');


%kdtree_delete(TextureCoords1_kdtree);
%kdtree_delete(TextureCoords2_kdtree);

if isfield(GM.Aux,'name') && isfield(GN.Aux,'name')
    rslt.Gname1 = GM.Aux.name;
    rslt.Gname2 = GN.Aux.name;
end
rslt.cPdist = cPdist;
rslt.cPmap = cPmap;
rslt.invcPmap = invcPmap;
rslt.TextureCoords1 = TextureCoords1;
rslt.TextureCoords2 = TextureCoords2;
rslt.translation = T;
rslt.orthogonal = R;

end
