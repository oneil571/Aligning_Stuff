%% Load excel file
try
    [~,~,raw] = xlsread(spreadsheetPath);
catch
    [~,~,raw] = xlsread(spreadsheetPath(1:end-1));
end
colsToDelete = [];
for i = 1:size(raw,2)
    if isnan(raw{1,i})
        colsToDelete = [colsToDelete i];
    end
end

raw(:,colsToDelete) = [];

workingPath = [projectDir specimenGroup '/'];

load([workingPath 'Names.mat'])



%% Map names to entries
nameMap = zeros(size(Names));
noMap = [];

for j = 2:size(raw,1)
    raw{j,1}(isspace(raw{j,1})) = [];
end
for i = 1:length(Names)
    for j = 2:size(raw,1)
        if contains(lower(Names{i}),lower(raw{j,1}))
            nameMap(i) = j;
            break
        elseif j == size(raw,1)
            noMap = [noMap i];
        end
    end
end

%% Extract all groupings of interest
Groups = containers.Map;

for g = 2:size(raw,2)
    curGroup = raw{1,g};
    curList = cell(size(Names));
    for i = 1:length(Names)
        if nameMap(i) == 0
            curList{i} = '';
            continue;
        end
        curList{i} = raw{nameMap(i),g};
        if isnan(curList{i})
            curList{i} = '';
        elseif isnumeric(curList{i})
            curList{i} = num2str(curList{i});
        end
    end
    Groups(curGroup) = curList;
end

save([workingPath 'Groups.mat'],'Groups');
save([workingPath 'groupsToCompare.mat'],'groupsToCompare');