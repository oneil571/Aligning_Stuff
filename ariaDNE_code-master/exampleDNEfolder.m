% This is an example computing ariaDNE for a user-specified folder




% set-up path:
clear; clc;

addpath('./utils');
addpath(genpath('./utils'));
%pathSetup('./alignedMeshes2');
addpath('./alignedMeshes2');
addpath(genpath('./alignedMeshes2'));

% pathSetup(BaseDirectory) %or provide a specified base directory


% compute ARIADNE
Options.distInfo = 'Geodeisic';
Options.cutThresh = 0;
foldername = ['./alignedMeshes2'];
bandwidth = 0.06;

result = computeDNEfolder(foldername, bandwidth, Options);
save('result.m', 'result')
fprintf('Computation is complete. Results are saved in result.m');

