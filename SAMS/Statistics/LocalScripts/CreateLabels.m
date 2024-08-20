if ~exist(infoPath)
    disp(['No data found, will ignore. Please ensure datapath is properly located or ' ...
      'please fill out labels by hand.']);
    return
end
[num,txt,stuff] = xlsread(infoPath);
load('D://Dropbox/TeethData/AnkleWrist/Ankles/Names.mat');
talusNames = Names;
talusGroups = cell(length(talusNames),size(txt,2)-2);
talusMap = containers.Map;
load('D://Dropbox/TeethData/AnkleWrist/Wrists/Names.mat');
hamateNames = Names;
hamateGroups = cell(length(hamateNames),size(txt,2)-2);
hamateMap = containers.Map;
for i = 2:100
    if strcmp(txt{i,2},'talus')
        ind = strFind(txt{i,1},talusNames);
        
        for j = 1:size(txt,2)-2
            talusGroups{ind,j} = txt{i,j+2};
        end
    else
        ind = strFind(txt{i,1},hamateNames);
        for j = 1:size(txt,2)-2
            hamateGroups{ind,j} = txt{i,j+2};
        end
    end
end
for i = 3:size(txt,2)
    hamateMap(txt{1,i}) = hamateGroups(:,i-2);
    talusMap(txt{1,i}) = talusGroups(:,i-2);
end
Groups = hamateMap;
save('D://Dropbox/TeethData/AnkleWrist/Wrists/Groups.mat','Groups');
Groups = talusMap;
save('D://Dropbox/TeethData/AnkleWrist/Ankles/Groups.mat','Groups');



