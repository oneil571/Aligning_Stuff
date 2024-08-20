function weightDist=ExtractMatchesPairsDistributions(temp,wtB,wtDecr,initMesh,...
    finalMesh,landmarks,neighborhoodSize,percDecr,minPerc,workingPath,...
    maxDepth)

% Script for computing the set of matched landmarks based off of forward
% propagation. Works by gradual matching, decreasing amount of sureness
% needed as time goes on in order to capture fuzzier correspondences. Extra
% means that curvature extrema are added to GP landmarks.
%
% Input:
% temp: bandwidth parameter for kernel.
% wtB: lower bound on minimum weight of path
% wtDecr: weight decrement in case wtB set too high
% initMesh/finalMesh: meshes to match landmarks for, arguments are
% symmetric
% landmarks: number of Gaussian process landmarks used in matching process.
% neighborhoodSize: size of neighborhood of GP landmark to count as match
% minPerc: initial minimum percent of probability mass required for match
% percDecr: decrement parameter for percent matching
% lowBndPerc: final percentage required.


%% Initialize variables, temporarily here until proper form figured out

%% Load necessary variables
    names_path = [workingPath 'Names.mat'];
    dist_path = [workingPath 'EOPExperimentData/dists_BD.mat'];
    maps_path = [workingPath 'EOPExperimentData/mapList_BD.mat'];
    flows_path = [workingPath 'EOPExperimentData/Flows_BD.mat'];
    samples_path = [workingPath 'ProcessedMAT/'];
%% Don't unnecessarily run
if initMesh == finalMesh
    error('Same mesh, nothing to do')
end

load(names_path);  
load(dist_path);

cPDistances = dists;

load(maps_path);
load(flows_path);

%% Declare and assign global variables to pass to compile results from internal recursion

%Relabel needs to be deprecated
global cPMap;           
cPMap = mapList;

%Number of landmarks to load/consider
global numLandmarks;    
numLandmarks = 18;
t =temp;         %temperature parameter for diffusion
       %matrix of cP distances
cPDistances = (cPDistances+cPDistances')/2;
meshes = cell(length(Names),1);
%% Load meshes
for i = 1:length(Names)
    if i == initMesh
        load([samples_path Names{i} '.mat']);
        meshes{i} = G;
        cPMap{i,i} = 1:G.nV;
    elseif i == finalMesh
        load([samples_path Names{i} '.mat']);
        meshes{i} = G;
        cPMap{i,i} = 1:G.nV;
    end
end
load([workingPath 'newLmkInds.mat']);

landmarksToTest_1 = newLmkInds{initMesh};
global vertWeight_12;      
vertWeight_12 = zeros(length(landmarksToTest_1),size(meshes{finalMesh}.V,2));


%% Load weighted flow matrices
curFlow_12 = Flows{initMesh,finalMesh};

weightedFlows_12 = sparse((cPDistances.^2)).*(curFlow_12);
for i = 1:size(weightedFlows_12,1)
    for j = 1:size(weightedFlows_12,2)
        if (weightedFlows_12(i,j) ~=0)
            weightedFlows_12(i,j) = exp(-weightedFlows_12(i,j)/t);
        end
    end
end
%% Run diffusion to gather points and distributions
while true
    DepthFirstSearchPlotting_12(weightedFlows_12,initMesh,finalMesh,landmarksToTest_1,1,wtB,1,initMesh,maxDepth);
    
    if norm(vertWeight_12)>0
        break;
    else %reset and redo computation with smaller weight
        vertWeight_12 = zeros(length(landmarksToTest_1),size(meshes{finalMesh}.V,2));
        wtB = wtB - wtDecr;
    end
end

weightDist = vertWeight_12;

end
function DepthFirstSearchPlotting_12(weightedFlows,init,final,lmks,weight,wtBnd,curDepth,curPath,maxDepth)
%global uniV2 convexComb totalWeight continueCounter
global cPMap
global vertWeight_12
global numLandmarks
    nextVerts = find(weightedFlows(init,:));
    
    for i = 1:length(nextVerts)
        if (curDepth == maxDepth) && (nextVerts(i) ~= final)
            continue;
        end
        nextLmks = cPMap{init,nextVerts(i)}(lmks);
        nextWeight = weight*weightedFlows(init,nextVerts(i));
        if (nextWeight < wtBnd) 
            continue;
        end
        
        nextPath = [curPath nextVerts(i)];
        if nextVerts(i) == final

            %scatter3(meshes{final}.V(1,nextLmks),meshes{final}.V(2,nextLmks),meshes{final}.V(3,nextLmks),40,repmat([1 1-(nextWeight-wtBnd)/(1-wtBnd) 1-(nextWeight-wtBnd)/(1-wtBnd)],length(nextLmks),1),'filled');
            for k = 1:numLandmarks
                vertWeight_12(k,nextLmks(k)) = vertWeight_12(k,nextLmks(k)) + nextWeight;
            end
        elseif curDepth < maxDepth
            DepthFirstSearchPlotting_12(weightedFlows,nextVerts(i),final,nextLmks,nextWeight,wtBnd,curDepth+1,nextPath,maxDepth);
        end
    end
            
end

