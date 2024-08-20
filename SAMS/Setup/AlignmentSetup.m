
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% setup parameters in this section 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Specify name of aligner
AlignerName = 'PuenteAlignment';

%% Specify input directory
meshesPath = '/gtmp/BoyerLab/trgao10/PNAS_HDM/';

%% Specify output directory
outputPath = '/gtmp/BoyerLab/trgao10/PNAS_HDM_output/';

%% Set parameters for the algorithm
restart = 1;                            %whether to flush previous output
iniNumPts = 128;                        %initial sampling for first pass
finNumPts = 256;                        %final sampling for second pass
ssType = 'FPS'; %%% 'FPS' | 'GPR'       %how to choose psuedolandmarks
                                        %FPS = farthest point sampling
                                        %GPR = Gaussian process landmarks
type = 'SDP'; %%% 'MST' | 'SPC' | 'SDP' %how to synchronize alignments
                                        %MST = minimum spanning tree
                                        %SPC = spectral relaxation
                                        %SDP = semidefinite program
allow_reflection = 1;                   % if true, allows for reflections
max_iter = 1000;                        % number of iterates for optimization

%% Set parameters for computational resources and notifications
use_cluster = 0;                        %Whether to use cluster; set
                                        %to 0 for use in personal machines
n_jobs = 120;                           %Number of jobs, not relevant for
                                        %personal use

%Email address to notify for completion in cluster job
email_notification = 'rravier@math.duke.edu'; 

