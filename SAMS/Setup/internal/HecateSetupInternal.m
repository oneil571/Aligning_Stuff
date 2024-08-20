%% DO NOT EDIT; SETS PATHS AND GUARANTEES PARAMETERS AS REQUIRED
workingPath = [projectDir specimenGroup '/'];
%Set and make relevant paths
warning('off','MATLAB:rmpath:DirNotFound');
rmpath(genpath([SAMSPath 'Alignment']));
rmpath(genpath([SAMSPath 'Mapping']));
path(path, genpath([SAMSPath 'Statistics/']));
path(path, genpath([SAMSPath 'utils/'])); 
path(path, genpath([SAMSPath 'VisualizationScripts']));
warning('on','MATLAB:rmpath:DirNotFound')
