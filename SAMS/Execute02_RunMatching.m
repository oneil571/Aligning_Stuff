initialize;
MappingSetup;
load([workingPath 'Flags.mat'])
load([workingPath 'Names.mat'])
addpath(genpath('./VisualizationScripts/'));
mappingPath = [workingPath 'MappingData/'];

MappingSetupAndFlowExtraction;

if ~isKey(Flags,'featureMappings') || ForceFeatureRecomputation
    ComputeFeatureMatching;
else
    disp('Feature mappings already computed');
    load([mappingPath 'FeatureMatches.mat']);
end

MatchesPairsOnFlyWrapper;
RefineInitialMatches;
disp(['Sparse correspondences computed. Please visualize with ' ...
    'visualizeLandmarkCorrespondences before continuing']);

featureList = cell(length(Names),1);
switch featureMap
    case 'Conf'
        for i = 1:length(Names)
            featureList{i} = meshList{i}.Aux.ConfMaxInds;
        end
    case 'Gauss'
        for i = 1:length(Names)
            featureList{i} = meshList{i}.Aux.GaussMaxInds;
        end
    case 'Mean'
        for i = 1:length(Names)
            featureList{i} = meshList{i}.Aux.MeanMinInds;
        end
    case 'DNE'
        for i = 1:length(Names)
            featureList{i} = meshList{i}.Aux.DNEMaxInds;
        end
end
