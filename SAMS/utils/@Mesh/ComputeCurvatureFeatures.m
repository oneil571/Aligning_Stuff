function Aux = ComputeCurvatureFeatures(GM,Aux,options)
%EXTRACTFEATURES Summary of this function goes here
%   Detailed explanation goes here

if nargin < 2
    options = struct();
end
SmoothCurvatureFields = getoptions(options,'SmoothCurvatureFields',3);
ConfMaxLocalWidth = getoptions(options,'ConfMaxLocalWidth',8);
GaussMaxLocalWidth = getoptions(options,'GaussMaxLocalWidth',10);
MeanMinLocalWidth = getoptions(options,'MeanMinLocalWidth',8);
DNEMaxLocalWidth = getoptions(options,'DNEMaxLocalWidth',8);
ExcludeBoundary = getoptions(options,'ExcludeBoundary',1);
Display = getoptions(options,'Display','off');

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% extract features (local maximum of conformal factors)
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[~,~,Cmin,Cmax,Cmean,Cgauss,~] = GM.ComputeCurvature(options);      %computes curvatures

%% Compute conformal factors
Lambda = CORR_calculate_conformal_factors(GM.F',GM.V',Aux.UniformizationV');
Aux.Conf = CORR_calculate_conformal_factors_verts(GM.F',Lambda);

nind = find(isnan(Aux.Conf)); %all NaN vertices
if(~isempty(nind))
    PP = GM.V'; %look for closest me-vert on original mesh which is not NaN later
    PP(nind,:)=[]; %remove the ones with NaN
    LL = Aux.Conf;
    LL(nind)=[];
    idxs = knnsearch(PP,GM.V(:,nind)');
    Aux.Conf(nind) = LL(idxs);
end

%% Smooth computed curvatures
Conf = Aux.Conf;
[~,TriAreas] = GM.ComputeSurfaceArea;

%Make weight matrix...
WeightMatrix = GM.F2V;
for j = 1:GM.nF
    WeightMatrix(j,:) = WeightMatrix(j,:).*(1./Aux.VertArea);
end
for j = 1:GM.nV
    WeightMatrix(:,j) = TriAreas.*WeightMatrix(:,j);
end

for j=1:SmoothCurvatureFields
    
    
    CgaussFace = mean(Cgauss(GM.F));
    Cgauss = CgaussFace*WeightMatrix;
    CmeanFace = mean(Cmean(GM.F));
    Cmean = CmeanFace*WeightMatrix;
    ConfFace = mean(Conf(GM.F));
    Conf = ConfFace*WeightMatrix;
end


Cgauss = Cgauss';
Aux.Conf = Conf';

DNE = Cmin.^2+Cmax.^2;
DNETruncInds = find(DNE>median(DNE));

[GaussMaxInds,~] = GM.FindLocalMax(Cgauss,GaussMaxLocalWidth,ExcludeBoundary);
[DNEMaxInds,~] = GM.FindLocalMax(Cmin.^2+Cmax.^2,DNEMaxLocalWidth,ExcludeBoundary);
[MeanMinInds,~] = GM.FindLocalMax(-abs(Cmean),MeanMinLocalWidth,ExcludeBoundary);
BB = GM.GetBoundingBox;
diam = sqrt(sum(abs(BB(:,2)-BB(:,1)).^2));
r = 0.3*diam;

Ring = GM.ComputeVertexRing;

if isfield(Aux,'Conf')
    [ConfMaxInds,~] = GM.FindLocalMax(Aux.Conf,ConfMaxLocalWidth,ExcludeBoundary);
    ToBeDelInds = [];
    for j=1:length(ConfMaxInds)
        if isempty(find(DNETruncInds == ConfMaxInds(j), 1))
            RingNBD = Ring{ConfMaxInds(j)};
            flag = 0;
            for k=1:length(RingNBD)
                if find(DNETruncInds == RingNBD(k))
                    flag = 1;
                    break;
                end
            end
            if (flag == 0)
                ToBeDelInds = [ToBeDelInds,j];
            end
        end
    end
    ConfMaxInds(ToBeDelInds) = [];
    
    Aux.ConfMaxInds = ConfMaxInds;
end

ToBeDelInds = [];
for j=1:length(GaussMaxInds)
    if isempty(find(DNETruncInds == GaussMaxInds(j), 1))
        ToBeDelInds = [ToBeDelInds,j];
    end
end
GaussMaxInds(ToBeDelInds) = [];

ToBeDelInds = [];
for j=1:length(DNEMaxInds)
    if isempty(find(DNETruncInds == DNEMaxInds(j), 1))
        ToBeDelInds = [ToBeDelInds,j];
    end
end
DNEMaxInds(ToBeDelInds) = [];

ToBeDelInds = [];
for j=1:length(MeanMinInds)
    if isempty(find(DNETruncInds == MeanMinInds(j), 1))
        ToBeDelInds = [ToBeDelInds,j];
    end
end
MeanMinInds(ToBeDelInds) = [];



% if strcmpi(Display, 'on')
%     if isfield(GM.Aux,'name')
%         figure('Name',['Features on ' GM.Aux.name]);
%     else
%         figure;
%     end
%     GM.draw();hold on;
%     set(gcf,'ToolBar','none');
%     if isfield(GM.Aux,'Conf')
%         scatter3(GM.V(1,ConfMaxInds),GM.V(2,ConfMaxInds),GM.V(3,ConfMaxInds),'g','filled');
%     end
%     scatter3(GM.V(1,GaussMaxInds),GM.V(2,GaussMaxInds),GM.V(3,GaussMaxInds),'r','filled');
%     scatter3(GM.V(1,GaussMinInds),GM.V(2,GaussMinInds),GM.V(3,GaussMinInds),'b','filled');
%     
%     set(gca, 'CameraUpVector', [0.8469,-0.5272,-0.0696]);
%     set(gca, 'CameraPosition', [0.0584,0.8255,-5.7263]);
%     set(gca, 'CameraTarget', [0.0122,-0.0075,0.0173]);
%     set(gca, 'CameraViewAngle', 10.5477);
% end

Aux.GaussMaxInds = GaussMaxInds;
Aux.MeanMinInds = MeanMinInds;
Aux.DNEMaxInds = DNEMaxInds;

end

