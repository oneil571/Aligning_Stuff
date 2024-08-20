
%% Set input and output paths, should change for every new dataset
%PLEASE MAKE SURE PATHS END IN SLASH
% Path of aligned input data. Below structure assumes PuenteAlignment output format
% OFF, PLY, OBJ
dataPath = '/Users/rileywilde/Auto3dgm_Python/alignedMeshes2';
projectDir = '/Users/rileywilde/SAMS/teeth_sep22';
addpath(dataPath)
addpath(genpath(dataPath))
numGPLmks = 200; %start with 200 for disc, 400 for sphere

%dataPath = 'C:\Users\rober\PersonalProjects\ToothAndClaw\ToothAndClawData\PNAS\teeth\aligned_output\aligned\';

distancePath = '';  %KEEP BLANK
%Auto3dgm: GPDists_High.mat
%Base path for everything in a project, may include multiple groups of
%nonhomolgous surfaces

%projectDir = 'C:\Users\rober\Dropbox\SAMSResults\SAMSDemo_ErrorTeeth\';
%The path of the current group you are working on, do not change if you are
%only working with one collection of homologous surfaces
specimenGroup = 'Default';
%Can ignore if only working with one set of homologous surfaces
%% Feature information. May need to tune.

%Force recomputation of features, default to 0 (false)
ForceFeatureRecomputation = 1;
%Number of Gaussian Process landmarks to draw; these are candidate
%pseudolandmarks. Choose based on collection, about 4-5 percent of number
%of vertices sufficient
   

%% Pseudolandmark matching parameters. Set and test based on visual output.

%Geometric features used to extract initial matches.
%Options are Conf = conformal factor, Gauss = Gaussian curvature, 
%Mean = absolute mean curvature, DNE = DNE.
%Gauss recommended to start, followed by Conf, DNE, Mean
featureMap = 'Gauss';
ForcePutativeMatching = 1;    %Force recomputation of putative matches
maxDepth = 3;           %Maximum number of links in path, keep between 3-5 for fast computations
maxNbrSize = 1;         %Amount of uncertainty tolerated in pseudolandmark matching
                        %Should be an integer between 1-3.


%% Parameters for final landmarking procedure
maxNumMatches = 12;     %Maximum number of landmarks matches used
                        %Set based on intuition for number of Type II/III
                        %landmarks usually required for analysis. Should be
                        %based on surface complexity.
minMatchDist = 0.18;     %Minimum distance needed between any landmarks to establish matching
                        %Landmarks not added if below this. Set negative if
                        %no error required.

%% Visualization Parameters
visNumRow = 8;          %When visualizing landmark correspondences,
                        %the number of surfaces to view at once.
%% Parameters for spherical Hecate only
numGPLmksHecate = 250;
featureMapHecate = 'Conf';
numLmksHecate = 15;
hecateFeatureInds = 10;

%% The following parameters are less likely to need to change

%Feature computation parameters, MUST BE INTEGER VALUED
%Amount of smoothing to apply to curvatures. Lower priority for tuning
SmoothCurvatureFields = 3;          %Amount of smoothing applied to curvatures          
ConfMaxLocalWidth = 8;              %Conformal factor maxima radius
GaussMaxLocalWidth = 10;            %Gauss curvature maxima radius
MeanMinLocalWidth = 8;              %Absolute mean curvature minima radius
DNEMaxLocalWidth = 8;               %DNE v1 maxima radius

%Putative Matching Parameters, 
minPerc = 0.5;          %Minimum relative likelihood needed to establish match.
                        %Set between 0.5 and 1
percDecr = 0.05;        %Decrement parameter for lowering required relative likelihood
                        %Should be positive, at most minPerc. Smalelr = more
                        %gradual decrement
pathWtDecr = .99;      %Decrement minimum weight of path considered
                        %Set between 0 and 1, numbers closer to one result
                        %in slower, higher precision computation
maxDistTol = 0.12;      %Possibly deprecated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DO NOT EDIT BELOW
MappingSetupInternal;

