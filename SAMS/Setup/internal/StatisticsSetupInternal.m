%% DO NOT EDIT; SETS PATHS AND GUARANTEES PARAMETERS AS REQUIRED

%Set and make relevant paths
warning('off','MATLAB:rmpath:DirNotFound');
rmpath(genpath([SAMSPath 'Alignment']));
rmpath(genpath([SAMSPath 'Mapping']));
path(path, genpath([SAMSPath 'Statistics/']));
path(path, genpath([SAMSPath 'utils/'])); 
path(path, genpath([SAMSPath 'VisualizationScripts']));
warning('on','MATLAB:rmpath:DirNotFound')
%% Start off by making directories for each group
clear k
metaDir = dir(projectDir);
SpecimenTypes = {};
if length(metaDir) == 3
    if exist([MetaGroupBasePath metaDir(3).name '/Groups.mat']) > 0
        SpecimenTypes = [SpecimenTypes metaDir(3).name];
        if ~exist('keys')
            load([MetaGroupBasePath metaDir(3).name '/Groups.mat']);
            keys = Groups.keys;
        end
    end
else
    for k = 3:length(metaDir)
        if isdir([projectDir metaDir(k).name])
            if exist([projectDir metaDir(k).name '/Groups.mat']) > 0
                SpecimenTypes = [SpecimenTypes metaDir(k).name];
                if ~exist('keys')
                    load([projectDir metaDir(k).name '/Groups.mat']);
                    keys = Groups.keys;
                end
            end


        end
    end
end
statPath = [projectDir 'Statistics/'];
interStatPath = [statPath 'InterGroup/'];
intraStatPath = [statPath 'IntraGroup/'];
touch(statPath); touch(interStatPath); touch(intraStatPath);
if exist('keys')
    for k = 1:length(keys)
        touch([interStatPath keys{k} '/']);
    end
end
touch([interStatPath 'Total/']);

TTestType = {'EqualCovariance','UnequalCovariance'};
if runPermutations
    TTestType = [TTestType {'EqualCovariance_Permutation','UnequalCovariance_Permutation'}];
end