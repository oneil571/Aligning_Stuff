load([workingPath 'GPDists.mat']);
load([workingPath 'Names.mat']);
load([workingPath 'MappingData/MatchesPairs_Thresheld.mat']);

rmpath(genpath([SAMSPath 'Mapping/external']));
addpath(genpath([SAMSPath 'utils/']));
%% Make directories if needed
orbDataPath = [workingPath 'OrbifoldData/'];
options.pointCloud = 0;
touch(orbDataPath);
meshList = cell(length(Names),1);
for i = 1:length(Names)
    load([workingPath 'ProcessedMAT/' Names{i} '.mat']);
    meshList{i} = G;
end
frechMean = find(min(sum(GPDists.^2))==sum(GPDists.^2));
for i = 1:length(Names)
    if i~=frechMean
        dirString = [orbDataPath Names{i} '__To__' Names{frechMean} '/'];
        if ~exist(dirString,'dir')
            mkdir(dirString);
        end
        meshList{i}.Write([dirString Names{i} '.off'],'off',options);
        meshList{frechMean}.Write([dirString Names{frechMean} '.off'],'off',options);
        fid = fopen([dirString Names{i} '.txt'],'w');
        frechid = fopen([dirString Names{frechMean} '.txt'],'w');
        curMatches = matchesPairs{i};
        for j = 1:size(curMatches,1)
            fprintf(fid,'%d\n',curMatches(j,1));
            fprintf(frechid,'%d\n',curMatches(j,2));
        end
        fclose(fid); fclose(frechid);
    end
end
