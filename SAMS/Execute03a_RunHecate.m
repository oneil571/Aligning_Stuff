
%% Initiailization
initialize;
HecateSetup;

%% Vertex
hecateDir = [workingPath 'HecateVertex/'];
touch(hecateDir)
clear H
CreateDiffusionMatrixVertex;

cfg.dirCollate = dirCollate;
cfg.colorSegments = colorSegments;
cfg.numMeshDisplay=1;

for numEigs = numEigsVec
    EigenDecomp;
    for numSegments = numSegmentsVec
        SpectralClustering;
        cfg.out = [hecateDir 'NumSegments_' num2str(numSegments) '/NumEigs_' num2str(numEigs) '/'];
        touch(cfg.out);
        segRes = SegResult(newMeshList, kIdx, vIdxCumSum);
        segRes.calc_data();
        segRes.export(cfg);
        close all;
    end
end

