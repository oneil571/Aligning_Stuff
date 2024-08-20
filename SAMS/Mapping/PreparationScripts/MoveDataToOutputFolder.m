%function starts filling in new working directory and extracts basic info
dataDir = dir(dataPath);
rawOFFPath = [workingPath '/RawOFF/'];
rawMATPath = [workingPath '/RawMAT/'];
badRawPath = [workingPath '/BadData/'];
touch(rawOFFPath);
touch(rawMATPath);
touch(badRawPath);
badOpenPath = [badRawPath 'ProblemsOpeningFile/'];
touch(badOpenPath);
Names = {}; options.pointCloud = 0;
progressbar
for i = 1:length(dataDir)
    progressbar((i-1)/(length(dataDir)))
    if length(dataDir(i).name) < 4
        continue;
    end
    try
        if strcmp(dataDir(i).name(end-3:end),'.off') ...
                || strcmp(dataDir(i).name(end-3:end),'.obj') ...
                || strcmp(dataDir(i).name(end-3:end),'.ply')
            Names = [Names dataDir(i).name(1:end-4)];
            switch dataDir(i).name(end-2:end)
                case 'off'
                    G = Mesh('off',[dataPath dataDir(i).name]);
                case 'obj'
                    G = Mesh('obj',[dataPath dataDir(i).name]);
                case 'ply'
                    G = Mesh('ply',[dataPath dataDir(i).name]);
                otherwise
                    continue;
            end
            G.Write([rawOFFPath dataDir(i).name(1:end-4) '.off'],'off',options);
            save([rawMATPath dataDir(i).name(1:end-4) '.mat'],'G');
        else
            error;
        end
    catch
        continue
    end
    
end
progressbar(1)
save([workingPath 'RawNames.mat'],'Names');

%This only works with Auto3dgm of old and should be deprecated

disp('NOTE: Loading prior distance matrix currently deprecated');
disp('Previous distances will not be used in final computations');
%if isfile(distancePath)
%    S = load(distancePath);
%    distName = fieldnames(S);
%    GPDists = getfield(S,distName{1});
%    Flags('hasDists') = 1;
%    save([workingPath 'GPDists.mat'],'GPDists');
%    save([workingPath 'Flags.mat'],'Flags');
%else
Flags('hasDists') = 0;
%end
save([workingPath 'Flags.mat'],'Flags');

