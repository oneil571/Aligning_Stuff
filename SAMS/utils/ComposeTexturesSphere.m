function [TextureCoords1,TextureCoords2] = ComposeTexturesSphere(TextureCoords1List,TextureCoords2List)
%COMPOSETEXTURECOORDSALONGPATH Summary of this function goes here
%   Detailed explanation goes here


TextureCoords1 = TextureCoords1List{end};
TextureCoords2 = TextureCoords2List{end};
for j = 2:length(TextureCoords1List)
    NextNodeTextureCoords1 = TextureCoords1List{end-j+1};
    NextNodeTextureCoords2 = TextureCoords2List{end-j+1};
    TextureCoords1 = PropagateTextureCoords(TextureCoords1,NextNodeTextureCoords1,NextNodeTextureCoords2);
end
end



function NewTextureCoords1 = PropagateTextureCoords(TextureCoords1,NextTextureCoords1,NextTextureCoords2)
% NextTextureCoords2 is the un-deformed version of TextureCoords1.
% This function returns NextTextureCoords1 under the same deformation that
% deformed NextTextureCoords2 to TextureCoords1.
% The returned NewTextureCoords1 should generically have different dimensions in
% column, compared to that of TextureCoords1.
%

BBox = [-1,-1,1,1;-1,1,-1,1];
NextTextureCoords2 = [NextTextureCoords2,BBox];
TextureCoords1 = [TextureCoords1,BBox];
Tri = delaunayTriangulation(NextTextureCoords2');

NewTextureCoords1 = zeros(size(NextTextureCoords1));
nV = size(NewTextureCoords1,2);
nF = size(Tri.ConnectivityList,1);
[new_ti,new_bc] = Tri.pointLocation(NextTextureCoords1');
for j = 1:nV
    
    if sum(isnan(new_bc(j,:))+isnan(new_bc(j,:))+isnan(new_bc(j,:))) == 0
        NewTextureCoords1(:,j) = TextureCoords1(:,Tri.ConnectivityList(new_ti(j),:))*new_bc(j,:)';
    else
        BC = Tri.cartesianToBarycentric((1:nF)',repmat(NextTextureCoords1(:,j)',nF,1));
        tind = find(all(BC>-1e-10,2));
        if(numel(tind)>1)
            tind = tind(1);
        elseif numel(tind)==0
            warning(['Point ', num2str(j), ' was not found in any triangle. It stays where it was.']);
            NewTextureCoords1(:,j) = NextTextureCoords1(:,j);
            continue;
        end
        BC = BC(tind,:);
        NewTextureCoords1(:,j) = TextureCoords1(:,Tri.ConnectivityList(tind,:))*BC';
    end
end

end

