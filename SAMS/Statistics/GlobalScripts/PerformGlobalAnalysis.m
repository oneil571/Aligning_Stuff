

%% Now do analysis
disp('Performing global analysis on MDS');
MDSGlobalAnalysis;
WeightedMeanMDSGlobalAnalysis;
disp('Performing global analysis on mean shapes');
SampleMeanGlobalAnalysis;
WeightedMeanShapeGlobalAnalysis;
disp('Performing global analysis on Frechet mean');
SampleFrechetMeanGlobalAnalysis;
WeightedFrechetMeanGlobalAnalysis;

disp('Writing all MDS Diagrams');
for i = 1:length(SpecimenTypes)
    curDir = [projectDir SpecimenTypes{i} '/'];
    PlotMDS(curDir,0); PlotMDS(curDir,1);
end
close all
