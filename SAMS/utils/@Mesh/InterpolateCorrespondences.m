function [rslt] = InterpolateCorrespondences(GM,GN,options)
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

%Coarse mapping: Nearest Neighbors

TPS_FEATURESM_INDS = options.GPMatches(:,1)';
TPS_FEATURESN_INDS = options.GPMatches(:,2)';


TextureCoords2 = GN.Aux.UniformizationV(1:2,:);
TextureCoords2(:,isnan(compl(TextureCoords2))) = ones(2,sum(isnan(compl(TextureCoords2))));
TextureCoords1 = GM.Aux.UniformizationV(1:2,:);
orig = TextureCoords1;

origFEATURESM = orig(:,TPS_FEATURESM_INDS);
baseFEATURESN = TextureCoords2(:,TPS_FEATURESN_INDS);
bestRot = eye(2);
bestDist = Inf;
for t = 0:pi/500:2*pi
    curRot = [cos(t) -sin(t);sin(t) cos(t)];
    curFeatures = origFEATURESM;
    for j = 1:size(origFEATURESM,2)
        curFeatures(:,j) = curRot*curFeatures(:,j);
    end
    if norm(curFeatures-baseFEATURESN,'fro') < bestDist
        bestRot = curRot;
        bestDist = norm(curFeatures-baseFEATURESN,'fro');
    end
end

for j = 1:size(TextureCoords1,2)
    TextureCoords1(:,j) = bestRot*TextureCoords1(:,j);
end
TPS_FEATURESM = DISCtoPLANE(TextureCoords1(:,TPS_FEATURESM_INDS)','d2p');
TPS_FEATURESN = DISCtoPLANE(TextureCoords2(:,TPS_FEATURESN_INDS)','d2p');

[TPS_FEATURESM,ia] = unique(TPS_FEATURESM,'rows');
TPS_FEATURESN = TPS_FEATURESN(ia,:);
[TPS_FEATURESN,ia] = unique(TPS_FEATURESN,'rows');
TPS_FEATURESM = TPS_FEATURESM(ia,:);

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
