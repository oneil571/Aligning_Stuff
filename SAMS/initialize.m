%% DO NOT EDIT
clear
SAMSPath = mfilename('fullpath');
SAMSPath = strsplit(SAMSPath,'initialize'); SAMSPath = SAMSPath{1};
path(pathdef);
path(path,genpath([SAMSPath 'Setup/']));
path(path,genpath([SAMSPath 'InternalScripts/']));