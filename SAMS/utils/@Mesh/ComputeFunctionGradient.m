function [vertGrad,faceGrad] = ComputeFunctionGradient(G,f)
%Computes the gradient of a function by definition in terms of directional
%derivative. Analogous to trigradient on MathWorks, but trigradient has
%some issues with matrix singularity, which shouldn't really occur...
%
%Input:
%G: Mesh
%f: function defined on vertices of a mesh
%
%Output:
%faceGrad: Gradient vectors on faces. Defined to be constant. 3xF matrix
%vertGrad: Gradient vectors on vertices. Defined to be weighted sum (in
%terms of area of neighboring faces) of gradients on faces. 3xV matrix.

%get face normals
[~,normalsFaces] = G.ComputeNormal;

%get face surface areas
[~,faceAreas]=G.ComputeSurfaceArea;

%check to make sure f is a row vector
if size(f,1)~=1
    f=f';
end

%make second set of vectors
secondVecs = repmat(f(G.F(1,:)),3,1).*(G.V(:,G.F(3,:))-G.V(:,G.F(2,:)))+...
    repmat(f(G.F(2,:)),3,1).*(G.V(:,G.F(1,:))-G.V(:,G.F(3,:)))+...
    repmat(f(G.F(3,:)),3,1).*(G.V(:,G.F(2,:))-G.V(:,G.F(1,:)));

%Compute face gradients
if size(faceAreas,1)~=1
    faceAreas = faceAreas';
end

faceGrad = cross(normalsFaces,secondVecs,1)./repmat(2*faceAreas,3,1);

vertFaceRing = CORR_compute_vertex_face_ring(G.F);
vertGrad = zeros(3,size(G.V,2));

for i = 1:size(G.V,2)
    vertGrad(:,i) = sum(repmat(faceAreas(vertFaceRing{i}),3,1)...
        .*faceGrad(:,vertFaceRing{i}),2)/sum(faceAreas(vertFaceRing{i}));
end
end

