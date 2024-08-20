

%% Set paths and organization of output, probably same as in MappingSetup
%Base path for everything in a project, may include multiple groups
projectDir = '/Users/rileywilde/SAMS/teeth_sep22/';
specimenGroup = 'Default';

                        
%% HDM/Hecate Settings
BNN = 5;                    %Number of nearest neighbors on base manifold
numEigsVec = [1:10];            %Number of eigenvectors to use in decompositions
                            %Highly variable data should use approximately
                            %5-10.
numSegmentsVec = [5];      %Vector of number of segments to consider
kMeansMaxIter = 5000;       %Number of iterations allowed for optimization to converge

%% Hecate Output Formatting
numMeshDisplay = 1;         %Number of specimens to display as example final segmentations
dirCollate = 0;             %Collate segment distribution groups, 0 = no, 1 = yes
colorSegments = 1;          %Whether to give each segment different colors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DO NOT EDIT BELOW
HecateSetupInternal;


