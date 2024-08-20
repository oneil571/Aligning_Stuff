function mData = isManifold( obj )
% ISMANIFOLD returns true if the mesh is a manifold, false otherwise.
%   input: mesh object
%   output:  a struct of the following fields:
% 
%       mData.manifold -  true if the mesh is a manifold,
%                         false otherwise.
% 
%       mData.nonManifoldEdges - a cell {i,j} of the indices of
%                                all the edges with more then 
%                                two adjacent faces.
% 
%       mData.nonManifoldVertices - a list of all the vertices
%                                   such that their 1-ring is 
%                                   not a disk.
% 
%       mData.nonManifoldBoundary - a list of all the boundary vertices
%                                   such that their 1-ring is 
%                                   not a disk. (aka they have more then
%                                   two neighbours on the boundary)
%
% Created by Nave Zelzer on march 31 2014.

% init step
Nf = obj.nF;
Nv = obj.nV;
F = obj.F;
F2V = obj.ComputeF2V;
boundary = obj.FindBoundaries;
mData.manifold = true;
mData.nonManifoldEdges = [];
mData.nonManifoldVertices = [];
mData.nonManifoldBoundary = [];

% we have three checks to do:
% 1 ) check that every edge only have one or two faces adjacent to it, no
%     more nor less!!!
% 
% 2 ) check that the mesh don't have an hourglass shape a.k.a exist a
%     vertex with 2-disk neighbourhood.
% 
% 3 ) check that the every boundary vertex have only two neighbours on the
%     boundary. aka the topology of the boundary is non disk shape.

% step 1:
% 1 ) check that every edge only have one or two faces adjacent to it, no
%     more nor less!!!
I = [F(1,:), F(2,:), F(3,:), F(2,:), F(3,:), F(1,:)];
J = [F(2,:), F(3,:), F(1,:), F(1,:), F(2,:), F(3,:)];
S = ones(1,6*Nf);
edges = sparse(I,J,S);
[i,j] = find(edges > 2);
if length(i) > 1
    mData.manifold = false;
    mData.nonManifoldEdges = {i,j};
end


% step 2:
% 2 ) check that the mesh don't have an hourglass shape a.k.a exist a
%     vertex with 2-disk neighbourhood.

for v=1:Nv
    %fprintf('v = %d\n',v);
    if ismember(v,boundary)         %Something is wrong with this method in this case
        continue;
    end
    first = 1;
    % find all adjacent faces to v
    adjf = obj.F(:,find(full(F2V(:,v))));
    % take the first face that's contain's exactly two boundary vertices
    % and get it's vertices
    b = find(sum(ismember(adjf,boundary)) == 2,1);
    % if we didn't find such face or we found one and v itself is not on
    % the boundary, then v is on a disk, else v is on a half disk.
    if isempty(b) || ~ismember(v,boundary)
        f = adjf(:,first);
        % take only the vertices on this face that are not v
        f = f(f ~= v);
        % take only the first one of those vertices
        f = f(1);   
    else   
        first = b;
        f = adjf(:,first);
        % take only the vertices on this face that are not v
        f = f(f ~= v);
        % take only the those vertices that are not on the boundary
        f = f(~ismember(f,boundary));
        % take only the first one of those vertices
        f = f(1);        
    end

    % find the only two possible faces that contain this vertex. aka the
    % face itself and the next neighbouring face.
    next = find(any(adjf==f));
    % take the next face
    next = next(next ~= first);
    %fprintf('%d->%d',first,next);
    counter = 1;
    while next ~= first 
        curr = next;
        currf = f;
        f = adjf(:,curr);
        f = f(f ~= v & f ~= currf);
        next = find(any(adjf==f));
        next = next(next ~= curr);
        if(isempty(next))
            next = first;
        else
            %fprintf('->%d',next);
        end
        counter = counter+1;
    end
    
    if(counter ~= size(adjf,2))
        mData.manifold = false;
        mData.nonManifoldVertices = [mData.nonManifoldVertices v];
    end
    %fprintf('\ncounter = %d, neighbours = %d\n\n',counter,size(adjf,2)); 
end

% step 3:
% 3 ) check that the every boundary vertex have only two neighbours on the
%     boundary. aka the topology of the boundary is non disk shape.
for b = boundary
    neigh = find(edges(b,:));
    neigh = find(ismember(neigh,boundary));
    if(length(neigh) > 2)
        mData.nonManifoldBoundary = [mData.nonManifoldBoundary b];
    end
end