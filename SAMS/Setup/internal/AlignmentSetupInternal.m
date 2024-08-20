%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DO NOT EDIT

%Find path of SAMS in system

%Clear path of previously existing directory
%Path of Aligner
SAMSPath= [SAMSPath 'Alignment/' AlignerName '/'];
rmpath([SAMSPath 'utils']); rmpath([SAMSPath 'Mapping']);
rmpath([SAMSPath 'Statistics']); rmpath([SAMSPath OrganizationScripts]);
path(pathdef);
path(path, genpath([SAMSPath 'software']), genpath([SAMSPath 'code']));

%Set mosek path
setenv('MOSEKLM_LICENSE_FILE', [SAMSPath 'software/mosek/mosek.lic'])