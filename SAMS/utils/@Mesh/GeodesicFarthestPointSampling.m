function [VertSampInd] = GeodesicFarthestPointSampling(G,SampNumb,InitialSamples)
% Geodesic farthest point sampling

if nargin<2
    SampNumb = G.nV;
end
if nargin<3
    rng('shuffle');
    InitialSamples = randi(G.nV,1,1);
else
    if size(InitialSamples,1)>size(InitialSamples,2)
        InitialSamples = InitialSamples';
    end
end

ProcessedSampNumb = length(InitialSamples);
VertSampInd = [InitialSamples, zeros(1,SampNumb-ProcessedSampNumb)];

%% Altered script to not use fast marching
meshAdj = sparse(pdist2(G.V',G.V').*G.A);
meshGraph = graph(meshAdj);
for k=(ProcessedSampNumb+1):SampNumb
    [~,VertSampInd(k)] = max(min(distances(meshGraph,VertSampInd(1:(k-1)))));
end
close all

end

