function [TextureCoords1,ImprMap] = TPSDeformationGP(GM,GN,InputMap,FeatureType,TextureCoords1,TextureCoords2,ref12,options)
%TPSDEFORMATION: TPS mesh deformation
%   Current design prefers geodesic mutually nearest neighboring features
%   than Euclidean mutually nearest neighboring features: if a feature is
%   assigned to different target features under these two ad-hoc
%   assignments, then the "Euclidean target" gets "wiped out" in the
%   "unique(TPS_FEATURESM_INDS)" section.

if nargin<7
    options = [];
end
%TextureCoords2_kdtree = getoptions(options,'TextureCoords2_kdtree',kdtree_build(TextureCoords2'));
GaussMinMatch = getoptions(options,'GaussMinMatch','on');

orig = TextureCoords1;

%%% extract matching features

TPS_FEATURESM_INDS = options.GPMatches(:,1);
TPS_FEATURESN_INDS = options.GPMatches(:,2);

%% Look at ConfMax to fill in unmatched details: best guess
[~,CONF_EUC_FEATURESN,preCONF_EUC_FEATURESM] = FindEuclideanMutuallyNearestNeighbors(GM,GN,InputMap,'ConfMax');
[~,CONF_FEATURESN,preCONF_FEATURESM] = FindMutuallyNearestNeighbors(GM,GN,InputMap,'ConfMax');
CONF_FEATURESM_INDS = [preCONF_FEATURESM;preCONF_EUC_FEATURESM];
CONF_FEATURESN_INDS = [CONF_FEATURESN;CONF_EUC_FEATURESN];
[CONF_FEATURESM_INDS,NoRepeatInds] = unique(CONF_FEATURESM_INDS);
CONF_FEATURESN_INDS = CONF_FEATURESN_INDS(NoRepeatInds);
% [~,TPS_GAUSSMIN_FEATURESN,preTPS_GAUSSMIN_FEATURESM] = FindMutuallyNearestNeighbors(GM,GN,InputMap,'GaussMin');
% CONF_FEATURESM_INDS = [CONF_FEATURESM_INDS;preTPS_GAUSSMIN_FEATURESM];
% CONF_FEATURESN_INDS = [CONF_FEATURESN_INDS;TPS_GAUSSMIN_FEATURESN];

%check to see if features are close to others already added; if so, don't
%add
[~,R,~] = MapToDist(GM.V,GN.V,InputMap,GM.Aux.VertArea);
tempTextureCoords1 = TextureCoords1;
tempTextureCoords2 = TextureCoords2;
if ref12 == 0
    tempTextureCoords1(2,:) = -tempTextureCoords1(2,:);
    tempTextureCoords2(2,:) = -tempTextureCoords2(2,:);
end

maxDist = .4;
matchedGPVerts_M = tempTextureCoords1(:,TPS_FEATURESM_INDS)'; matchedGPVerts_N = tempTextureCoords2(:,TPS_FEATURESN_INDS)';
[~,dists_M] = knnsearch(matchedGPVerts_M,tempTextureCoords1(:,CONF_FEATURESM_INDS)');
[~,dists_N] = knnsearch(matchedGPVerts_N,tempTextureCoords2(:,CONF_FEATURESN_INDS)');

% for i = 1:length(CONF_FEATURESM_INDS)
%     if (dists_M(i) > .5*maxDist) && (dists_N(i) > .5*maxDist) && ...
%             (norm(tempTextureCoords1(:,CONF_FEATURESM_INDS(i))-...
%             tempTextureCoords2(:,CONF_FEATURESN_INDS(i))) < maxDist)
%             TPS_FEATURESM_INDS o= [TPS_FEATURESM_INDS;CONF_FEATURESM_INDS(i)];
%             TPS_FEATURESN_INDS = [TPS_FEATURESN_INDS;CONF_FEATURESN_INDS(i)];
%     end
% end

%% Add GP features
TPS_FEATURESM = DISCtoPLANE(TextureCoords1(:,TPS_FEATURESM_INDS)','d2p');
TPS_FEATURESN = DISCtoPLANE(TextureCoords2(:,TPS_FEATURESN_INDS)','d2p');

[TPS_FEATURESM,ia] = unique(TPS_FEATURESM,'rows');
TPS_FEATURESN = TPS_FEATURESN(ia,:);
[TPS_FEATURESN,ia] = unique(TPS_FEATURESN,'rows');
TPS_FEATURESM = TPS_FEATURESM(ia,:);





%% Interpolate
if (length(TPS_FEATURESM)>3) % TPS (Thin Plate Spline)
    tP = DISCtoPLANE(TextureCoords1','d2p');
    [ftps] = TEETH_calc_tps(TPS_FEATURESM,TPS_FEATURESN-TPS_FEATURESM);
    pt = tP + TEETH_eval_tps(ftps,tP);
    TextureCoords1 = DISCtoPLANE(pt,'p2d')';
%     TextureCoords1(:,isnan(compl(TextureCoords1))) = 1;
elseif (length(TPS_FEATURESM)==3) % Affine Transformation
    tP = DISCtoPLANE(TextureCoords1','d2p');
    [A,b] = PlanarThreePtsDeform(TPS_FEATURESM,TPS_FEATURESN);
    pt = [A,b]*[tP';ones(1,size(tP,1))];
    TextureCoords1 = DISCtoPLANE(pt','p2d')';
%     TextureCoords1(:,isnan(compl(TextureCoords1))) = 1;
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

%%% construct vertex permutation map
ImprMap = knnsearch(TextureCoords2', TextureCoords1');

end

